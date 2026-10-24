"""Tests for Unified Training & Cloud Execution Telemetry Logger and GCP Account Policy Guard (Job 442)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from factorforge.telemetry.training_logger import (
    REQUIRED_GCP_ACCOUNT,
    GCPAccountPolicyError,
    TrainingTelemetryLogger,
    verify_gcp_account_policy,
)


def test_gcp_account_policy_blocks_legacy_free_credit_project_and_unauthorized_accounts() -> None:
    """Ensures unauthorized accounts and blocked cdx projects are rejected."""
    with pytest.raises(GCPAccountPolicyError, match="not authorized"):
        verify_gcp_account_policy(
            project_id="fair-conduit-510507-p5",
            account_email="cdx-user@gmail.com",
        )

    with pytest.raises(GCPAccountPolicyError, match="blocked"):
        verify_gcp_account_policy(
            project_id="cdx-legacy-free-tier",
            account_email=REQUIRED_GCP_ACCOUNT,
        )

    with pytest.raises(GCPAccountPolicyError, match="bucket"):
        verify_gcp_account_policy(
            project_id="fair-conduit-510507-p5",
            account_email=REQUIRED_GCP_ACCOUNT,
            bucket_name="cdx-legacy-bucket",
        )


def test_gcp_account_policy_accepts_authorized_munkyu_eijex_configuration(
    tmp_path: Path,
) -> None:
    """Verifies that munkyu.kim@eijex.com (fair-conduit-510507-p5) passes policy checks."""
    key_file = tmp_path / "factorforge-runner-key.json"
    key_file.write_text(
        json.dumps(
            {
                "type": "service_account",
                "project_id": "fair-conduit-510507-p5",
                "client_email": "factorforge-runner@fair-conduit-510507-p5.iam.gserviceaccount.com",
            }
        ),
        encoding="utf-8",
    )

    verified = verify_gcp_account_policy(
        project_id="fair-conduit-510507-p5",
        account_email="munkyu.kim@eijex.com",
        key_path=key_file,
        bucket_name="factorforge-fair-conduit",
    )
    assert verified["account_email"] == "munkyu.kim@eijex.com"
    assert verified["project_id"] == "fair-conduit-510507-p5"
    assert verified["bucket_name"] == "factorforge-fair-conduit"


def test_training_telemetry_logger_persists_hash_only_records_and_sqlite_tables(
    tmp_path: Path,
) -> None:
    """Verifies JSONL and SQLite persistence across 1DP, 3DP, 4DP, local_rewriter, and sllm_training."""
    logger = TrainingTelemetryLogger(tmp_path / "telemetry")
    protein = "MKELVA"
    dna = "ATGAAAGAATTGGTGGCG"

    rec_1dp = logger.log_optimization_run(
        engine_mode="1DP",
        protein_sequence=protein,
        coding_dna=dna,
        score=0.033804,
        gc_percent=44.44,
        cai_global=0.8125,
        is_feasible=True,
        step_evaluations=120,
        wall_time_ms=1.5,
        subsystem="factorforge",
    )
    rec_4dp = logger.log_optimization_run(
        engine_mode="4DP",
        protein_sequence=protein,
        coding_dna=dna,
        score=0.0215,
        gc_percent=44.44,
        cai_global=0.795,
        is_feasible=True,
        step_evaluations=120,
        wall_time_ms=1.8,
        extra_metrics={"tunnel_length_assumption": 35, "half_width": 2},
        subsystem="agentos",
    )
    rec_rw = logger.log_optimization_run(
        engine_mode="local_rewriter",
        protein_sequence=protein,
        coding_dna=dna,
        score=0.033804,
        gc_percent=44.44,
        cai_global=0.8125,
        is_feasible=True,
        step_evaluations=45,
        wall_time_ms=0.9,
        extra_metrics={"accepted_rewrites": 4, "optimality_gap_vs_dp": 0.0},
        subsystem="validationhub",
    )
    rec_ml = logger.log_ml_training_stage(
        stage_id="008-sllm-balanced-training-review",
        corpus_digest="a" * 64,
        train_sequences=11,
        val_sequences=3,
        test_sequences=3,
        val_nll=1.1816,
        test_nll=1.1243,
        baseline_test_nll=1.1169,
    )

    assert len(rec_1dp.canonical_sequence_sha256) == 64
    assert rec_4dp.evidence_tier == "UNVERIFIED_PUBLIC_CORPUS"
    assert rec_rw.biological_validation == "NOT_AVAILABLE"
    assert rec_ml.gcp_target_account == "munkyu.kim@eijex.com"

    # Confirm no raw sequence is written anywhere in the JSONL file
    raw_jsonl = logger.jsonl_path.read_text(encoding="utf-8")
    assert dna not in raw_jsonl
    assert protein not in raw_jsonl

    counts = logger.summary_counts()
    assert counts["execution_runs_count"] == 4
    assert counts["computational_assessments_count"] == 3
    assert counts["training_telemetry_count"] == 4
    assert counts["by_engine_mode"] == {
        "1DP": 1,
        "4DP": 1,
        "local_rewriter": 1,
        "sllm_training": 1,
    }


def test_training_telemetry_logger_rejects_raw_sequence_leak_or_translation_mismatch(
    tmp_path: Path,
) -> None:
    """Ensures raw sequences in extra_metrics or non-synonymous DNA are blocked fail-closed."""
    logger = TrainingTelemetryLogger(tmp_path / "telemetry")

    with pytest.raises(ValueError, match="does not translate identically"):
        logger.log_optimization_run(
            engine_mode="1DP",
            protein_sequence="MKELVA",
            coding_dna="ATGTTTTTTTTTTTTTTT",
            score=0.0,
            gc_percent=20.0,
            cai_global=0.5,
            is_feasible=True,
        )

    with pytest.raises(ValueError, match="Raw nucleotide sequence|Prohibited raw sequence"):
        logger.log_optimization_run(
            engine_mode="1DP",
            protein_sequence="MKELVA",
            coding_dna="ATGAAAGAATTGGTGGCG",
            score=0.0,
            gc_percent=44.44,
            cai_global=0.8,
            is_feasible=True,
            extra_metrics={"leaked_dna": "ATGAAAGAATTGGTGGCG"},
        )


@pytest.mark.parametrize(
    "extra",
    [
        {"score": 0},
        {"is_feasible": False},
        {"note": "private observation"},
        {"value": float("nan")},
        {"nested": {"x": 1}},
    ],
)
def test_extra_metrics_integrity(tmp_path, extra):
    logger = TrainingTelemetryLogger(tmp_path)
    with pytest.raises(ValueError):
        logger.log_optimization_run(
            engine_mode="1DP",
            protein_sequence="MKELVA",
            coding_dna="ATGAAAGAATTGGTGGCG",
            score=0,
            gc_percent=44.44,
            cai_global=0.8,
            is_feasible=True,
            extra_metrics=extra,
        )
    assert logger.summary_counts()["training_telemetry_count"] == 0
    assert not logger.jsonl_path.exists()


@pytest.mark.parametrize(
    "overrides",
    [
        {"score": float("inf")},
        {"is_feasible": "false"},
        {"step_evaluations": 1.5},
        {"wall_time_ms": -1},
    ],
)
def test_invalid_core_metrics(tmp_path, overrides):
    kwargs = {
        "engine_mode": "1DP",
        "protein_sequence": "MKELVA",
        "coding_dna": "ATGAAAGAATTGGTGGCG",
        "score": 0,
        "gc_percent": 44.44,
        "cai_global": 0.8,
        "is_feasible": True,
    }
    kwargs.update(overrides)
    with pytest.raises(ValueError):
        TrainingTelemetryLogger(tmp_path).log_optimization_run(**kwargs)


def test_missing_credentials_are_not_verified(tmp_path):
    with pytest.raises(GCPAccountPolicyError, match="unavailable"):
        verify_gcp_account_policy(key_path=tmp_path / "missing.json")


def test_ml_digest_and_baseline_integrity(tmp_path):
    logger = TrainingTelemetryLogger(tmp_path)
    kwargs = {
        "stage_id": "008-review",
        "corpus_digest": "a" * 64,
        "train_sequences": 1,
        "val_sequences": 1,
        "test_sequences": 1,
        "val_nll": 1.0,
        "test_nll": 1.2,
        "baseline_test_nll": 1.0,
    }
    with pytest.raises(ValueError, match="cannot override"):
        logger.log_ml_training_stage(**kwargs, extra_metrics={"beats_baseline": True})
    kwargs["corpus_digest"] = "label-is-not-a-content-hash"
    with pytest.raises(ValueError, match="SHA-256"):
        logger.log_ml_training_stage(**kwargs)


@pytest.mark.parametrize("field", ["engine_id", "engine_version"])
def test_sequence_like_metadata_is_rejected(tmp_path, field):
    kwargs = {
        "engine_mode": "1DP",
        "protein_sequence": "MKELVA",
        "coding_dna": "ATGAAAGAATTGGTGGCG",
        "score": 0,
        "gc_percent": 44.44,
        "cai_global": 0.8,
        "is_feasible": True,
    }
    kwargs[field] = "MELKVMELKVMELKVMELKV"
    logger = TrainingTelemetryLogger(tmp_path)
    with pytest.raises(ValueError, match="protein-like"):
        logger.log_optimization_run(**kwargs)
    assert logger.summary_counts()["training_telemetry_count"] == 0
