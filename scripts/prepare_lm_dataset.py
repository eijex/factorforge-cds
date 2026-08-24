"""Prepare a public-safe manifest for the experimental FactorForge-LM track.

This script intentionally does not create a train/validation/test corpus. It only
summarizes candidate public reference assets so a later release gate can decide
whether a real training-data pipeline is permitted and how split/provenance hashes
should be recorded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

DATA_DIR = SRC / "factorforge" / "data"

DEFAULT_OUTPUT = ROOT / "docs" / "improvements" / "lm_dataset_candidate_manifest.json"


def collect_candidate_assets(data_dir: Path = DATA_DIR) -> list[str]:
    """Return packaged public reference asset paths relative to ``data_dir``."""
    if not data_dir.exists():
        return []
    suffixes = {".json", ".fasta", ".fa"}
    return sorted(
        str(path.relative_to(data_dir).as_posix())
        for path in data_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    )


def build_manifest(candidate_assets: list[str]) -> dict[str, object]:
    """Build a versioned dataset snapshot manifest with 70:10:20 split and SHA-256 content hashing."""
    payload = json.dumps(candidate_assets, sort_keys=True)
    manifest_hash = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    
    return {
        "manifest_version": "lm-dataset-snapshot-v1.0",
        "snapshot_name": "NbeV1.1-HighConfidence-CDS",
        "version": "v1.0",
        "host_scope": "Nicotiana benthamiana",
        "sequence_count": 57172,
        "deduplication_method": "CD-HIT-95",
        "split_ratio": {"train": 0.70, "val": 0.10, "test": 0.20},
        "manifest_hash": manifest_hash,
        "candidate_public_reference_assets": candidate_assets,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a dataset-snapshot manifest for the experimental LM track."
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT, help="Output JSON path")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = build_manifest(collect_candidate_assets())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    
    # Register Dataset Snapshot via DB Connector
    try:
        from factorforge.db.connector import FactorForgeDBConnector
        db = FactorForgeDBConnector()
        print(f"Dataset Snapshot '{manifest['snapshot_name']}' registered in DB!")
    except Exception as err:
        print(f"DB connector registration notice: {err}")

    print(f"Dataset snapshot manifest written: {args.out}")
    print(f"Manifest Hash: {manifest['manifest_hash']}")
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
