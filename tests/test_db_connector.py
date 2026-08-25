"""PyTest suite for the local Context-Aware FactorForgeDBConnector."""

import sqlite3
import tempfile
from contextlib import closing
from pathlib import Path

import pytest

from factorforge.db.connector import FactorForgeDBConnector


def test_db_connector_save_candidate_and_evaluations():
    with tempfile.TemporaryDirectory() as tmp:
        db_path = str(Path(tmp) / "factorforge.db")

        db = FactorForgeDBConnector(db_path=db_path)
        res = db.save_candidate_with_evaluations(
            candidate_data={
                "engine": "lm",
                "model_version": "v3.5.0-SynCodonLM-V2",
                "optimized_sequence": "ATGGCTAAATGGTAA",
                "cai": 0.84,
                "gc_percent": 41.2,
                "type2is_clean": True,
            },
            evaluations=[
                {"constraint_code": "AA_IDENTITY", "status": "PASS", "observed_value": 1.0},
                {"constraint_code": "TYPE_IIS_CLEAN", "status": "PASS", "observed_value": 0},
            ],
            construct_id="PUBLIC-SYNTHETIC-TEST",
        )

        assert res["package_id"] >= 1
        assert res["candidate_id"] >= 1
        assert len(res["frozen_manifest_hash"]) == 64


def _candidate_asset_manifest(db: FactorForgeDBConnector) -> dict[str, object]:
    manifest: dict[str, object] = {
        "manifest_version": "lm-candidate-assets-v1",
        "snapshot_name": "FactorForge-Packaged-Public-Reference-Assets",
        "version": "v1.0",
        "status": "candidate-assets-only-no-training-split",
        "host_scope": "mixed-packaged-reference-assets",
        "sequence_count": None,
        "deduplication_method": None,
        "split_ratio": None,
        "candidate_public_reference_assets": ["example.json"],
        "asset_fingerprints": [
            {"path": "example.json", "size_bytes": 3, "sha256": "a" * 64}
        ],
    }
    manifest["manifest_hash"] = db._manifest_hash(manifest)
    return manifest


def test_register_dataset_snapshot_is_real_and_idempotent():
    with tempfile.TemporaryDirectory() as tmp:
        db_path = str(Path(tmp) / "factorforge.db")
        db = FactorForgeDBConnector(db_path=db_path)
        manifest = _candidate_asset_manifest(db)

        first_id = db.register_dataset_snapshot(manifest)
        second_id = db.register_dataset_snapshot(manifest)

        assert first_id == second_id
        with closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute(
                "SELECT COUNT(*), status, sequence_count "
                "FROM factorforge_dataset_snapshots"
            ).fetchone()
        assert row == (1, "candidate-assets-only-no-training-split", None)


def test_register_dataset_snapshot_rejects_hash_mismatch():
    with tempfile.TemporaryDirectory() as tmp:
        db = FactorForgeDBConnector(db_path=str(Path(tmp) / "factorforge.db"))
        manifest = _candidate_asset_manifest(db)
        manifest["manifest_hash"] = "0" * 64

        with pytest.raises(ValueError, match="hash mismatch"):
            db.register_dataset_snapshot(manifest)
