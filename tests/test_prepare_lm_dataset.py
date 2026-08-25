"""Regression tests for the public-safe LM candidate-asset manifest."""

import json
import sqlite3
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_prepare_lm_dataset_records_only_verified_candidate_assets(tmp_path):
    output_path = tmp_path / "candidate-assets.json"
    db_path = tmp_path / "factorforge.db"

    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "prepare_lm_dataset.py"),
            "--out",
            str(output_path),
            "--db-path",
            str(db_path),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "candidate-assets-only-no-training-split"
    assert manifest["sequence_count"] is None
    assert manifest["deduplication_method"] is None
    assert manifest["split_ratio"] is None
    assert manifest["asset_fingerprints"]
    assert all(len(item["sha256"]) == 64 for item in manifest["asset_fingerprints"])
    assert "registered in local SQLite" in completed.stdout

    with sqlite3.connect(db_path) as conn:
        stored = conn.execute(
            "SELECT snapshot_name, status, manifest_hash "
            "FROM factorforge_dataset_snapshots"
        ).fetchone()
    assert stored == (
        manifest["snapshot_name"],
        manifest["status"],
        manifest["manifest_hash"],
    )
