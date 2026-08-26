"""FactorForge Database Connector Module (v0.3 Postgres PostgreSQL Track).

Implements a Context-Aware, Policy-Driven Rule Engine connector that connects to the
canonical v0.3 Tri-Tier architecture:
1. common: Canonical Sequence Registry & Biological Context
2. factorforge: Rule Engine & Execution Provenance
3. validationhub: Physical Evidence & Claims
4. agentops: Policy Enforcement
5. ml: Research Scaffolds
"""

import hashlib
import json
import os
import sqlite3
import uuid
from contextlib import closing
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import psycopg2
from psycopg2.extras import DictCursor

from factorforge.utils.sequence_identity import canonicalize_sequence


_DB_PATH_UNSET = object()


class FactorForgeDBConnector:
    """Python DB connector for FactorForge Context-Aware Rule Engine."""

    def __init__(
        self,
        dsn: Optional[str] = None,
        db_path: Optional[str] | object = _DB_PATH_UNSET,
    ):
        """Initialize a PostgreSQL connector or an explicit local SQLite checkpoint.

        ``db_path`` is retained for candidate-asset research checkpoints and
        tests. Omitting it selects PostgreSQL; explicitly passing ``db_path``
        (including ``None``) selects the local SQLite checkpoint contract.
        """
        self._sqlite_mode = db_path is not _DB_PATH_UNSET
        if self._sqlite_mode:
            if db_path is None:
                db_dir = Path(__file__).resolve().parents[3] / "data" / "db"
                db_dir.mkdir(parents=True, exist_ok=True)
                self.db_path = str(db_dir / "factorforge_relational.db")
            else:
                self.db_path = str(db_path)
            self._init_sqlite_checkpoint()
            return
        self.dsn = dsn or os.environ.get(
            "FACTORFORGE_DB_DSN", 
            "host=127.0.0.1 port=5432 user=postgres password=postgres dbname=postgres"
        )
        self._init_db_extensions()

    def _sqlite_connection(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        return conn

    def _init_sqlite_checkpoint(self) -> None:
        """Create the minimal normalized schema for local research checkpoints."""
        with closing(self._sqlite_connection()) as conn:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS factorforge_sequence_registry (
                    sequence_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sequence_hash TEXT UNIQUE NOT NULL,
                    molecule_type TEXT NOT NULL,
                    sequence_length INTEGER NOT NULL,
                    raw_sequence TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                );
                CREATE TABLE IF NOT EXISTS factorforge_candidates (
                    candidate_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sequence_id INTEGER NOT NULL REFERENCES factorforge_sequence_registry(sequence_id),
                    engine_name TEXT NOT NULL,
                    model_version TEXT NOT NULL,
                    overall_status TEXT NOT NULL,
                    cai_value REAL,
                    gc_percent REAL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                );
                CREATE TABLE IF NOT EXISTS factorforge_constraint_evaluations (
                    evaluation_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    candidate_id INTEGER NOT NULL REFERENCES factorforge_candidates(candidate_id),
                    constraint_code TEXT NOT NULL,
                    status TEXT NOT NULL,
                    observed_value REAL,
                    details_json TEXT
                );
                CREATE TABLE IF NOT EXISTS factorforge_design_packages (
                    package_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    candidate_id INTEGER NOT NULL REFERENCES factorforge_candidates(candidate_id),
                    construct_id TEXT NOT NULL,
                    frozen_manifest_hash TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                );
                CREATE TABLE IF NOT EXISTS factorforge_dataset_snapshots (
                    snapshot_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    snapshot_name TEXT NOT NULL,
                    version TEXT NOT NULL,
                    status TEXT NOT NULL,
                    host_scope TEXT,
                    sequence_count INTEGER,
                    deduplication_method TEXT,
                    split_ratio_json TEXT,
                    manifest_hash TEXT NOT NULL,
                    manifest_json TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    UNIQUE(snapshot_name, version, manifest_hash)
                );
                """
            )
            conn.commit()

    @contextmanager
    def get_connection(self):
        """Context manager for DB connections."""
        conn = psycopg2.connect(self.dsn)
        try:
            yield conn
        finally:
            conn.close()

    def _init_db_extensions(self):
        """Ensure necessary PostgreSQL extensions/structures exist if needed."""
        # Our DDL (schema.sql) is now the authority. We do not CREATE TABLE here.
        # This connector merely operates on the existing v0.3 schema.
        pass

    @staticmethod
    def _manifest_hash(manifest: Mapping[str, Any]) -> str:
        """Return the canonical SHA-256 hash for a manifest payload."""
        payload = dict(manifest)
        payload.pop("manifest_hash", None)
        canonical = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def register_dataset_snapshot(self, manifest: Mapping[str, Any]) -> str:
        """Register verified manifest metadata in the ml.dataset_snapshots research checkpoint."""
        required = {"snapshot_name", "version", "status", "manifest_hash"}
        missing = sorted(required.difference(manifest))
        if missing:
            raise ValueError(f"Dataset snapshot manifest missing fields: {missing}")

        expected_hash = self._manifest_hash(manifest)
        supplied_hash = str(manifest["manifest_hash"])
        if supplied_hash != expected_hash:
            raise ValueError(
                "Dataset snapshot manifest hash mismatch: "
                f"expected {expected_hash}, received {supplied_hash}"
            )

        if self._sqlite_mode:
            manifest_json = json.dumps(
                dict(manifest), sort_keys=True, separators=(",", ":"), ensure_ascii=False
            )
            split_ratio = manifest.get("split_ratio")
            split_ratio_json = (
                json.dumps(split_ratio, sort_keys=True, separators=(",", ":"))
                if split_ratio is not None
                else None
            )
            with closing(self._sqlite_connection()) as conn:
                conn.execute(
                    """
                    INSERT OR IGNORE INTO factorforge_dataset_snapshots
                    (snapshot_name, version, status, host_scope, sequence_count,
                     deduplication_method, split_ratio_json, manifest_hash, manifest_json)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        manifest["snapshot_name"],
                        manifest["version"],
                        manifest["status"],
                        manifest.get("host_scope"),
                        manifest.get("sequence_count"),
                        manifest.get("deduplication_method"),
                        split_ratio_json,
                        supplied_hash,
                        manifest_json,
                    ),
                )
                row = conn.execute(
                    """
                    SELECT snapshot_id FROM factorforge_dataset_snapshots
                    WHERE snapshot_name = ? AND version = ? AND manifest_hash = ?
                    """,
                    (manifest["snapshot_name"], manifest["version"], supplied_hash),
                ).fetchone()
                if row is None:
                    raise RuntimeError("Dataset snapshot registration did not return a row")
                conn.commit()
                return str(row["snapshot_id"])

        snapshot_id = str(uuid.uuid4())
        
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                # Check for existing
                cursor.execute(
                    """
                    SELECT snapshot_id
                    FROM ml.dataset_snapshots
                    WHERE snapshot_name = %s AND version = %s AND manifest_hash = %s
                    """,
                    (manifest["snapshot_name"], manifest["version"], supplied_hash),
                )
                row = cursor.fetchone()
                if row:
                    return str(row[0])

                cursor.execute(
                    """
                    INSERT INTO ml.dataset_snapshots
                    (snapshot_id, snapshot_name, version, status, manifest_hash)
                    VALUES (%s, %s, %s, %s, %s)
                    """,
                    (
                        snapshot_id,
                        manifest["snapshot_name"],
                        manifest["version"],
                        manifest["status"],
                        supplied_hash,
                    ),
                )
            conn.commit()
            return snapshot_id

    def register_sequence(self, raw_sequence: str, molecule_class: str = "CDS") -> str:
        """Registers a canonical sequence in common.sequences.
        
        Note: Requires an underlying common.artifacts entry as per v0.3 spec invariant.
        """
        normalized_seq, seq_hash = canonicalize_sequence(molecule_class, raw_sequence)

        if self._sqlite_mode:
            with closing(self._sqlite_connection()) as conn:
                row = conn.execute(
                    "SELECT sequence_id FROM factorforge_sequence_registry WHERE sequence_hash = ?",
                    (seq_hash,),
                ).fetchone()
                if row:
                    return str(row["sequence_id"])
                cursor = conn.execute(
                    """
                    INSERT INTO factorforge_sequence_registry
                    (sequence_hash, molecule_type, sequence_length, raw_sequence)
                    VALUES (?, ?, ?, ?)
                    """,
                    (seq_hash, molecule_class, len(normalized_seq), normalized_seq),
                )
                conn.commit()
                return str(cursor.lastrowid)
        
        with self.get_connection() as conn:
            with conn.cursor(cursor_factory=DictCursor) as cursor:
                cursor.execute(
                    "SELECT sequence_id FROM common.sequences WHERE canonical_sequence_sha256 = %s",
                    (seq_hash,),
                )
                row = cursor.fetchone()
                if row:
                    return str(row["sequence_id"])

                # Insert dummy artifact to satisfy FK invariant (in real app, this points to actual fasta file in GCS/S3)
                artifact_id = str(uuid.uuid4())
                cursor.execute(
                    """
                    INSERT INTO common.artifacts (artifact_id, artifact_uri, sha256_hash, artifact_type)
                    VALUES (%s, %s, %s, 'fasta')
                    """,
                    (artifact_id, f"internal://generated/{seq_hash}.fasta", seq_hash)
                )

                seq_id = str(uuid.uuid4())
                cursor.execute(
                    """
                    INSERT INTO common.sequences (sequence_id, molecule_class, canonical_sequence_sha256, artifact_id, length)
                    VALUES (%s, %s, %s, %s, %s)
                    """,
                    (seq_id, molecule_class, seq_hash, artifact_id, len(normalized_seq)),
                )
            conn.commit()
            return seq_id

    def save_candidate_with_evaluations(
        self,
        candidate_data: Dict[str, Any],
        evaluations: List[Dict[str, Any]],
        construct_id: str = "CF-CANDIDATE-001",
    ) -> Dict[str, Any]:
        """Persist a candidate in the explicit local SQLite checkpoint."""
        if not self._sqlite_mode:
            raise RuntimeError(
                "save_candidate_with_evaluations is a local SQLite checkpoint operation"
            )
        sequence_id = int(self.register_sequence(candidate_data["optimized_sequence"]))
        with closing(self._sqlite_connection()) as conn:
            cursor = conn.execute(
                """
                INSERT INTO factorforge_candidates
                (sequence_id, engine_name, model_version, overall_status, cai_value, gc_percent)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    sequence_id,
                    candidate_data.get("engine", "lm"),
                    candidate_data.get("model_version", "v3.5.0-SynCodonLM-V2"),
                    "PASS" if candidate_data.get("type2is_clean", True) else "REJECT",
                    candidate_data.get("cai"),
                    candidate_data.get("gc_percent"),
                ),
            )
            candidate_id = int(cursor.lastrowid)
            for evaluation in evaluations:
                conn.execute(
                    """
                    INSERT INTO factorforge_constraint_evaluations
                    (candidate_id, constraint_code, status, observed_value, details_json)
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    (
                        candidate_id,
                        evaluation["constraint_code"],
                        evaluation["status"],
                        evaluation.get("observed_value"),
                        json.dumps(evaluation.get("details", {})),
                    ),
                )
            manifest_payload = json.dumps(
                {
                    "candidate_id": candidate_id,
                    "sequence_id": sequence_id,
                    "construct_id": construct_id,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            manifest_hash = hashlib.sha256(manifest_payload.encode("utf-8")).hexdigest()
            package_cursor = conn.execute(
                """
                INSERT INTO factorforge_design_packages
                (candidate_id, construct_id, frozen_manifest_hash)
                VALUES (?, ?, ?)
                """,
                (candidate_id, construct_id, manifest_hash),
            )
            conn.commit()
            return {
                "package_id": int(package_cursor.lastrowid),
                "candidate_id": candidate_id,
                "sequence_id": sequence_id,
                "frozen_manifest_hash": manifest_hash,
            }

    def save_computational_provenance(
        self,
        run_metadata: Dict[str, Any],
        candidate_data: Dict[str, Any],
        check_results: List[Dict[str, Any]],
    ) -> str:
        """
        Atomically saves the execution provenance (Phase 1):
        1. design_run
        2. candidate
        3. check_results
        
        Rolls back the entire transaction if any step fails. Does NOT create a design_package.
        """
        raw_seq = candidate_data["optimized_sequence"]
        # Sequence registration could technically be part of the transaction, but it's idempotent.
        seq_id = self.register_sequence(raw_seq)
        
        run_id = str(uuid.uuid4())
        candidate_id = str(uuid.uuid4())

        # The context manager automatically commits on success, and rolls back on exception
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                # 1. Insert Design Run
                cursor.execute(
                    """
                    INSERT INTO factorforge.design_runs 
                    (design_run_id, execution_origin, actual_engine_name, actual_profile_name, 
                     generation_performed, analysis_mode, evidence_level, runner_entrypoint)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                    """,
                    (
                        run_id,
                        run_metadata.get("execution_origin", "standalone_cli"),
                        run_metadata["actual_engine_name"],
                        run_metadata.get("actual_profile_name"),
                        run_metadata.get("generation_performed", True),
                        run_metadata.get("analysis_mode", "cli_generation"),
                        run_metadata.get("evidence_level", "prospective_factorforge"),
                        run_metadata.get("runner_entrypoint", "factorforge.cli.main"),
                    )
                )

                # 2. Insert Candidate (Metrics)
                overall_status = candidate_data.get("computational_status", "passed")
                cursor.execute(
                    """
                    INSERT INTO factorforge.candidates
                    (candidate_id, design_run_id, sequence_id, computational_status, cai_value, gc_percent)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    """,
                    (
                        candidate_id,
                        run_id,
                        seq_id,
                        overall_status,
                        candidate_data.get("cai", None),
                        candidate_data.get("gc_percent", None),
                    ),
                )

                # 3. Insert Check Results (Constraints)
                for chk in check_results:
                    check_result_id = str(uuid.uuid4())
                    cursor.execute(
                        """
                        INSERT INTO factorforge.check_results
                        (check_result_id, candidate_id, check_type, result, observed_value, details_json)
                        VALUES (%s, %s, %s, %s, %s, %s)
                        """,
                        (
                            check_result_id,
                            candidate_id,
                            chk["check_domain"],  # e.g., sequence_integrity, advisory_sequence_risk
                            chk["result"],        # PASS, WARN, FAIL
                            chk.get("observed_value"),
                            json.dumps(chk.get("details_json", {})),
                        ),
                    )
            # Automatic commit here; if exception occurs, it rolls back.
            conn.commit()

        return run_id

    def save_ml_evaluation_provenance(
        self,
        evaluation_run_data: Dict[str, Any],
        evaluation_checks: List[Dict[str, Any]],
    ) -> str:
        """
        Atomically saves the ML Evaluation Integrity provenance:
        1. evaluation_runs
        2. evaluation_checks
        """
        run_id = str(uuid.uuid4())
        
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                # 1. Insert Evaluation Run
                cursor.execute(
                    """
                    INSERT INTO ml.evaluation_runs
                    (evaluation_run_id, model_id, evaluation_snapshot_id, engine_name, evaluation_protocol, status)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    """,
                    (
                        run_id,
                        evaluation_run_data.get("model_id"),
                        evaluation_run_data.get("evaluation_snapshot_id"),
                        evaluation_run_data.get("engine_name", "UNKNOWN"),
                        evaluation_run_data.get("evaluation_protocol", "standard_benchmark"),
                        evaluation_run_data.get("status", "completed"),
                    )
                )

                # 2. Insert Evaluation Checks
                for chk in evaluation_checks:
                    check_id = str(uuid.uuid4())
                    cursor.execute(
                        """
                        INSERT INTO ml.evaluation_checks
                        (check_id, evaluation_run_id, check_type, result, threshold, match_count, details_json)
                        VALUES (%s, %s, %s, %s, %s, %s, %s)
                        """,
                        (
                            check_id,
                            run_id,
                            chk["check_type"],
                            chk["result"],
                            chk.get("threshold"),
                            chk.get("match_count"),
                            json.dumps(chk.get("details_json", {})),
                        )
                    )
            conn.commit()
            
        return run_id
