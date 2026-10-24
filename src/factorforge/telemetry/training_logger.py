"""Unified Training & Cloud Execution Telemetry Logger (Job 442).

Bridges FactorForge (1DP / 3DP / 4DP / ConstrainedLocalRewriter / sLLM) runs with
AgentOS execution traces, ResearchOS sequence assessments, and ValidationHub
computational evidence under strict Hash-Only Sequence Governance and GCP
account isolation (`munkyu.kim@eijex.com`).

Governance Invariants:
1. GCP Account & Project Guard:
   Checks the configured project, bucket and service-account metadata before
   runner client creation. This is not proof of an authenticated human identity.
2. Hash-Only Sequence Storage:
   Never writes raw DNA/mRNA/protein sequences into telemetry JSONL or SQLite tables.
   Only normalized 64-character SHA-256 digests (`canonical_sequence_sha256`,
   `aa_sha256`) and numerical features are persisted.
3. Truthful Evidence Metadata:
   Every computational telemetry record explicitly carries
   `evidence_tier="UNVERIFIED_PUBLIC_CORPUS"` and
   `biological_validation="NOT_AVAILABLE"`.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import sqlite3
import uuid
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from factorforge.analysis.metrics import translate_dna

REQUIRED_GCP_ACCOUNT = "munkyu.kim@eijex.com"
AUTHORIZED_GCP_PROJECT_ID = "fair-conduit-510507-p5"
AUTHORIZED_GCS_BUCKET = "factorforge-fair-conduit"
AUTHORIZED_KEY_PATH = str(Path.home() / ".gcp" / "factorforge-runner-key.json")
BLOCKED_LEGACY_GCP_PROJECTS = frozenset({"cdx-legacy-free-tier", "unauthorized-free-tier"})
BLOCKED_LEGACY_GCS_BUCKETS = frozenset({"cdx-legacy-bucket"})

_DNA_SEQ_RE = re.compile(r"^[ACGTU]{15,}$", re.IGNORECASE)
_PROTEIN_SEQ_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{15,}$", re.IGNORECASE)
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class GCPAccountPolicyError(RuntimeError):
    """Raised when a GCP operation targets a blocked project or lacks munkyu.kim@eijex.com credentials."""


def verify_gcp_account_policy(
    project_id: str | None = None,
    account_email: str | None = None,
    key_path: str | Path | None = None,
    bucket_name: str | None = None,
) -> dict[str, str]:
    """Verifies that GCP configuration targets the munkyu.kim@eijex.com workspace (fair-conduit-510507-p5)."""
    resolved_account = (
        account_email or os.environ.get("EIJEX_GCP_ACCOUNT") or REQUIRED_GCP_ACCOUNT
    ).strip()
    if resolved_account.lower() != REQUIRED_GCP_ACCOUNT:
        raise GCPAccountPolicyError(
            f"GCP account '{resolved_account}' is not authorized. "
            f"All GCP training and storage must use '{REQUIRED_GCP_ACCOUNT}'."
        )

    resolved_project = (
        project_id or os.environ.get("EIJEX_GCP_PROJECT_ID") or AUTHORIZED_GCP_PROJECT_ID
    ).strip()
    if (
        resolved_project.lower() in BLOCKED_LEGACY_GCP_PROJECTS
        or resolved_project != AUTHORIZED_GCP_PROJECT_ID
    ):
        raise GCPAccountPolicyError(
            f"GCP project '{resolved_project}' is blocked (depleted free-credit account). "
            f"Use '{AUTHORIZED_GCP_PROJECT_ID}' under '{REQUIRED_GCP_ACCOUNT}'."
        )

    resolved_bucket = (
        bucket_name or os.environ.get("EIJEX_GCP_BUCKET") or AUTHORIZED_GCS_BUCKET
    ).strip()
    if resolved_bucket != AUTHORIZED_GCS_BUCKET:
        raise GCPAccountPolicyError(
            f"GCS bucket '{resolved_bucket}' belongs to a blocked project and cannot be used."
        )

    resolved_key = (
        str(key_path)
        if key_path is not None
        else os.environ.get("EIJEX_GCP_KEY_PATH", AUTHORIZED_KEY_PATH)
    ).strip()

    key_file = Path(resolved_key)
    if not key_file.is_file():
        raise GCPAccountPolicyError("Configured service-account credential file is unavailable")
    if key_file.is_file():
        try:
            key_data = json.loads(key_file.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise GCPAccountPolicyError("Failed to parse configured credential JSON") from exc
        if not isinstance(key_data, dict):
            raise GCPAccountPolicyError("Credential metadata must be a JSON object")
        expected_email = f"factorforge-runner@{resolved_project}.iam.gserviceaccount.com"
        if (
            key_data.get("type") != "service_account"
            or key_data.get("client_email") != expected_email
        ):
            raise GCPAccountPolicyError(
                "Credential service-account metadata does not match configured runner"
            )
        file_proj = str(key_data.get("project_id", "")).strip()
        if file_proj.lower() in BLOCKED_LEGACY_GCP_PROJECTS:
            raise GCPAccountPolicyError(
                f"Credential file {key_file} points to blocked project '{file_proj}'."
            )
        if file_proj != resolved_project:
            raise GCPAccountPolicyError(
                f"Credential project_id '{file_proj}' does not match target project '{resolved_project}'."
            )

    return {
        "account_email": resolved_account,
        "project_id": resolved_project,
        "bucket_name": resolved_bucket,
        "key_path": str(key_file),
    }


def _sha256_hex(text: str) -> str:
    norm = "".join(text.upper().split())
    return hashlib.sha256(norm.encode("utf-8")).hexdigest()


def _merge_metrics(base: dict[str, Any], extra: dict[str, Any] | None) -> None:
    """Accept numeric observations only; callers cannot replace authoritative fields."""
    if extra is None:
        return
    if not isinstance(extra, dict):
        raise TypeError("extra_metrics must be a dictionary")
    _assert_no_raw_sequences(extra, "extra_metrics")
    for key, value in extra.items():
        if not isinstance(key, str) or not re.fullmatch(r"[a-z_][a-z0-9_]{0,63}", key):
            raise ValueError("extra_metrics keys must be lowercase metric identifiers")
        if key in base:
            raise ValueError("extra_metrics cannot override authoritative metrics")
        if not isinstance(value, (int, float, bool)) or not math.isfinite(value):
            raise ValueError("extra_metrics values must be finite numeric observations")
    base.update(extra)


def _finite_number(value: Any, name: str, minimum: float | None = None) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")


def _assert_no_raw_sequences(payload: Any, path: str = "root") -> None:
    """Recursively verifies that no raw nucleotide/protein sequences leak into telemetry records."""
    if isinstance(payload, dict):
        for k, v in payload.items():
            _assert_no_raw_sequences(str(k), "field_name")
            key_lower = str(k).lower()
            if key_lower in {
                "sequence",
                "coding_dna",
                "raw_sequence",
                "protein_sequence",
                "initial_feasible_sequence",
                "exhaustive_best_sequence",
                "dp_sequence",
                "rewriter_sequence",
            }:
                raise ValueError(
                    f"Prohibited raw sequence field '{path}.{k}' in telemetry payload."
                )
            _assert_no_raw_sequences(v, f"{path}.{k}")
    elif isinstance(payload, (list, tuple)):
        for idx, item in enumerate(payload):
            _assert_no_raw_sequences(item, f"{path}[{idx}]")
    elif isinstance(payload, str) and not _SHA256_RE.fullmatch(payload):
        normalized = "".join(payload.split())
        if _DNA_SEQ_RE.fullmatch(normalized) or _PROTEIN_SEQ_RE.fullmatch(normalized):
            raise ValueError(
                f"Raw nucleotide sequence or protein-like text detected at '{path}' in telemetry payload."
            )


@dataclass(frozen=True)
class TelemetryRecord:
    """Hash-only structured telemetry record for training & orchestration audit."""

    run_id: str
    task_id: str
    timestamp_utc: str
    subsystem: str  # "factorforge" | "agentos" | "researchos" | "validationhub"
    engine_mode: str  # "1DP" | "3DP" | "4DP" | "local_rewriter" | "sllm_training"
    engine_id: str
    engine_version: str
    subject_digest: str
    canonical_sequence_sha256: str
    aa_sha256: str
    protein_length: int
    total_nt: int
    metrics: dict[str, Any]
    gcp_target_account: str
    gcp_sync_status: str
    evidence_tier: str
    biological_validation: str


class TrainingTelemetryLogger:
    """Append-only JSONL and SQLite telemetry logger bridging FactorForge, AgentOS, ResearchOS, and ValidationHub."""

    def __init__(self, telemetry_dir: Path | str) -> None:
        self.telemetry_dir = Path(telemetry_dir).resolve()
        self.telemetry_dir.mkdir(parents=True, exist_ok=True)
        self.jsonl_path = self.telemetry_dir / "training_execution_ledger.jsonl"
        self.sqlite_path = self.telemetry_dir / "execution_trace.sqlite"
        self._init_sqlite()

    @contextmanager
    def _connection(self):
        db = sqlite3.connect(self.sqlite_path, timeout=30)
        db.row_factory = sqlite3.Row
        try:
            db.execute("PRAGMA synchronous=FULL")
            with db:
                yield db
        finally:
            db.close()

    def _init_sqlite(self) -> None:
        with self._connection() as db:
            # 1. AgentOS-compatible execution_runs table
            db.execute(
                """CREATE TABLE IF NOT EXISTS execution_runs (
                    run_id TEXT PRIMARY KEY,
                    task_id TEXT UNIQUE NOT NULL,
                    subject_digest TEXT NOT NULL,
                    status TEXT NOT NULL,
                    started_at TEXT NOT NULL,
                    finished_at TEXT,
                    error_type TEXT
                )"""
            )
            # 2. ValidationHub-compatible computational_assessments table
            db.execute(
                """CREATE TABLE IF NOT EXISTS computational_assessments (
                    task_id TEXT PRIMARY KEY,
                    assessment_id TEXT UNIQUE NOT NULL,
                    run_id TEXT NOT NULL,
                    subject_digest TEXT NOT NULL,
                    schema_version TEXT NOT NULL,
                    evidence_json TEXT NOT NULL
                )"""
            )
            # 3. Unified ML & DP training feature table (hash-only sequences)
            db.execute(
                """CREATE TABLE IF NOT EXISTS training_telemetry (
                    run_id TEXT PRIMARY KEY,
                    task_id TEXT UNIQUE NOT NULL,
                    timestamp_utc TEXT NOT NULL,
                    subsystem TEXT NOT NULL,
                    engine_mode TEXT NOT NULL,
                    engine_id TEXT NOT NULL,
                    engine_version TEXT NOT NULL,
                    subject_digest TEXT NOT NULL,
                    canonical_sequence_sha256 TEXT NOT NULL,
                    aa_sha256 TEXT NOT NULL,
                    protein_length INTEGER NOT NULL,
                    total_nt INTEGER NOT NULL,
                    metrics_json TEXT NOT NULL,
                    gcp_target_account TEXT NOT NULL,
                    gcp_sync_status TEXT NOT NULL,
                    evidence_tier TEXT NOT NULL,
                    biological_validation TEXT NOT NULL
                )"""
            )

    def log_optimization_run(
        self,
        *,
        engine_mode: str,
        protein_sequence: str,
        coding_dna: str,
        score: float,
        gc_percent: float,
        cai_global: float,
        is_feasible: bool,
        step_evaluations: int = 0,
        wall_time_ms: float = 0.0,
        extra_metrics: dict[str, Any] | None = None,
        subsystem: str = "factorforge",
        engine_id: str = "dp_v2_3_0",
        engine_version: str = "2.3.0",
        task_id: str | None = None,
    ) -> TelemetryRecord:
        """Logs a 1DP / 3DP / 4DP / LocalRewriter optimization run in hash-only format."""
        if subsystem not in {"factorforge", "agentos", "researchos", "validationhub"}:
            raise ValueError(f"Unsupported subsystem '{subsystem}'")
        if engine_mode not in {"1DP", "3DP", "4DP", "local_rewriter", "sllm_training"}:
            raise ValueError(f"Unsupported engine_mode '{engine_mode}'")

        _finite_number(score, "score")
        _finite_number(gc_percent, "gc_percent", 0)
        _finite_number(cai_global, "cai_global", 0)
        _finite_number(wall_time_ms, "wall_time_ms", 0)
        if gc_percent > 100 or cai_global > 1:
            raise ValueError("gc_percent or cai_global outside metric bounds")
        if type(is_feasible) is not bool:
            raise ValueError("is_feasible must be a boolean")
        if type(step_evaluations) is not int or step_evaluations < 0:
            raise ValueError("step_evaluations must be a nonnegative integer")
        clean_prot = "".join(protein_sequence.upper().split()).removesuffix("*")
        clean_dna = "".join(coding_dna.upper().split())
        if (
            not clean_prot
            or len(clean_dna) != len(clean_prot) * 3
            or translate_dna(clean_dna) != clean_prot
        ):
            raise ValueError(
                "Telemetry rejected: coding_dna does not translate identically to protein_sequence."
            )

        seq_sha256 = _sha256_hex(clean_dna)
        prot_sha256 = _sha256_hex(clean_prot)

        now_iso = datetime.now(timezone.utc).isoformat()
        resolved_task_id = str(uuid.UUID(str(task_id))) if task_id else str(uuid.uuid4())
        run_id = str(uuid.uuid4())

        metrics_payload: dict[str, Any] = {
            "score": round(float(score), 6),
            "gc_percent": round(float(gc_percent), 4),
            "cai_global": round(float(cai_global), 4),
            "is_feasible": bool(is_feasible),
            "step_evaluations": int(step_evaluations),
            "wall_time_ms": round(float(wall_time_ms), 3),
        }
        _merge_metrics(metrics_payload, extra_metrics)

        subject_obj = {
            "subsystem": subsystem,
            "engine_mode": engine_mode,
            "engine_id": engine_id,
            "engine_version": engine_version,
            "aa_sha256": prot_sha256,
            "canonical_sequence_sha256": seq_sha256,
        }
        subject_digest = hashlib.sha256(
            json.dumps(subject_obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()

        record = TelemetryRecord(
            run_id=run_id,
            task_id=resolved_task_id,
            timestamp_utc=now_iso,
            subsystem=subsystem,
            engine_mode=engine_mode,
            engine_id=engine_id,
            engine_version=engine_version,
            subject_digest=subject_digest,
            canonical_sequence_sha256=seq_sha256,
            aa_sha256=prot_sha256,
            protein_length=len(clean_prot),
            total_nt=len(clean_dna),
            metrics=metrics_payload,
            gcp_target_account=REQUIRED_GCP_ACCOUNT,
            gcp_sync_status="PENDING_GCP_MUNKYU_ACCOUNT_SYNC",
            evidence_tier="UNVERIFIED_PUBLIC_CORPUS",
            biological_validation="NOT_AVAILABLE",
        )

        rec_dict = asdict(record)
        _assert_no_raw_sequences(rec_dict)

        assessment_evidence = {
            "analysis_mode": (
                "baseline_design" if engine_mode == "1DP" else "structure_conditioned_design"
            ),
            "engine_mode": engine_mode,
            "generation_performed": True,
            "validation_passed": bool(is_feasible),
            "cai_global": round(float(cai_global), 4),
            "gc_percent": round(float(gc_percent), 4),
            "score": round(float(score), 6),
            "step_evaluations": int(step_evaluations),
            "canonical_sequence_sha256": seq_sha256,
            "aa_sha256": prot_sha256,
            "evidence_tier": "UNVERIFIED_PUBLIC_CORPUS",
            "biological_validation": "NOT_AVAILABLE",
        }

        with self._connection() as db:
            db.execute(
                "INSERT INTO execution_runs VALUES (?, ?, ?, ?, ?, ?, NULL)",
                (
                    run_id,
                    resolved_task_id,
                    subject_digest,
                    "succeeded" if is_feasible else "failed",
                    now_iso,
                    now_iso,
                ),
            )
            db.execute(
                "INSERT INTO computational_assessments VALUES (?, ?, ?, ?, ?, ?)",
                (
                    resolved_task_id,
                    str(uuid.uuid4()),
                    run_id,
                    subject_digest,
                    "0.1.0",
                    json.dumps(assessment_evidence, sort_keys=True),
                ),
            )
            db.execute(
                """INSERT INTO training_telemetry VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )""",
                (
                    record.run_id,
                    record.task_id,
                    record.timestamp_utc,
                    record.subsystem,
                    record.engine_mode,
                    record.engine_id,
                    record.engine_version,
                    record.subject_digest,
                    record.canonical_sequence_sha256,
                    record.aa_sha256,
                    record.protein_length,
                    record.total_nt,
                    json.dumps(record.metrics, sort_keys=True),
                    record.gcp_target_account,
                    record.gcp_sync_status,
                    record.evidence_tier,
                    record.biological_validation,
                ),
            )

        with open(self.jsonl_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec_dict, ensure_ascii=False, sort_keys=True) + "\n")

        return record

    def log_ml_training_stage(
        self,
        *,
        stage_id: str,
        corpus_digest: str,
        train_sequences: int,
        val_sequences: int,
        test_sequences: int,
        val_nll: float,
        test_nll: float,
        baseline_test_nll: float,
        extra_metrics: dict[str, Any] | None = None,
        task_id: str | None = None,
    ) -> TelemetryRecord:
        """Logs an sLLM / ML training stage summary for longitudinal tracking."""
        if not isinstance(corpus_digest, str) or not _SHA256_RE.fullmatch(corpus_digest):
            raise ValueError("corpus_digest must be a SHA-256 digest, not a corpus label")
        for name, value in (
            ("train_sequences", train_sequences),
            ("val_sequences", val_sequences),
            ("test_sequences", test_sequences),
        ):
            if type(value) is not int or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        for name, value in (
            ("val_nll", val_nll),
            ("test_nll", test_nll),
            ("baseline_test_nll", baseline_test_nll),
        ):
            _finite_number(value, name, 0)

        now_iso = datetime.now(timezone.utc).isoformat()
        resolved_task_id = str(uuid.UUID(str(task_id))) if task_id else str(uuid.uuid4())
        run_id = str(uuid.uuid4())

        metrics_payload: dict[str, Any] = {
            "stage_id": stage_id,
            "train_sequences": int(train_sequences),
            "val_sequences": int(val_sequences),
            "test_sequences": int(test_sequences),
            "val_nll": round(float(val_nll), 6),
            "test_nll": round(float(test_nll), 6),
            "baseline_test_nll": round(float(baseline_test_nll), 6),
            "beats_baseline": bool(test_nll < baseline_test_nll),
        }
        _merge_metrics(metrics_payload, extra_metrics)

        subject_obj = {
            "subsystem": "researchos",
            "engine_mode": "sllm_training",
            "stage_id": stage_id,
            "corpus_digest": corpus_digest,
        }
        subject_digest = hashlib.sha256(
            json.dumps(subject_obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()

        record = TelemetryRecord(
            run_id=run_id,
            task_id=resolved_task_id,
            timestamp_utc=now_iso,
            subsystem="researchos",
            engine_mode="sllm_training",
            engine_id="sllm_codon_lm",
            engine_version=stage_id,
            subject_digest=subject_digest,
            canonical_sequence_sha256=corpus_digest,
            aa_sha256=corpus_digest,
            protein_length=0,
            total_nt=0,
            metrics=metrics_payload,
            gcp_target_account=REQUIRED_GCP_ACCOUNT,
            gcp_sync_status="PENDING_GCP_MUNKYU_ACCOUNT_SYNC",
            evidence_tier="UNVERIFIED_PUBLIC_CORPUS",
            biological_validation="NOT_AVAILABLE",
        )

        rec_dict = asdict(record)
        _assert_no_raw_sequences(rec_dict)

        with self._connection() as db:
            db.execute(
                "INSERT INTO execution_runs VALUES (?, ?, ?, ?, ?, ?, NULL)",
                (run_id, resolved_task_id, subject_digest, "succeeded", now_iso, now_iso),
            )
            db.execute(
                """INSERT INTO training_telemetry VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )""",
                (
                    record.run_id,
                    record.task_id,
                    record.timestamp_utc,
                    record.subsystem,
                    record.engine_mode,
                    record.engine_id,
                    record.engine_version,
                    record.subject_digest,
                    record.canonical_sequence_sha256,
                    record.aa_sha256,
                    record.protein_length,
                    record.total_nt,
                    json.dumps(record.metrics, sort_keys=True),
                    record.gcp_target_account,
                    record.gcp_sync_status,
                    record.evidence_tier,
                    record.biological_validation,
                ),
            )

        with open(self.jsonl_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec_dict, ensure_ascii=False, sort_keys=True) + "\n")

        return record

    def list_records(self, engine_mode: str | None = None) -> list[dict[str, Any]]:
        """Returns all logged telemetry records, optionally filtered by engine_mode."""
        if not self.jsonl_path.exists():
            return []
        rows: list[dict[str, Any]] = []
        with open(self.jsonl_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if engine_mode is None or rec.get("engine_mode") == engine_mode:
                    rows.append(rec)
        return rows

    def summary_counts(self) -> dict[str, Any]:
        """Returns aggregate counts across SQLite and JSONL stores."""
        with self._connection() as db:
            exec_count = db.execute("SELECT COUNT(*) FROM execution_runs").fetchone()[0]
            assess_count = db.execute("SELECT COUNT(*) FROM computational_assessments").fetchone()[
                0
            ]
            telem_count = db.execute("SELECT COUNT(*) FROM training_telemetry").fetchone()[0]
            by_mode_rows = db.execute(
                "SELECT engine_mode, COUNT(*) AS cnt FROM training_telemetry GROUP BY engine_mode"
            ).fetchall()
        by_mode = {r["engine_mode"]: r["cnt"] for r in by_mode_rows}
        return {
            "execution_runs_count": exec_count,
            "computational_assessments_count": assess_count,
            "training_telemetry_count": telem_count,
            "by_engine_mode": by_mode,
            "gcp_target_account": REQUIRED_GCP_ACCOUNT,
        }
