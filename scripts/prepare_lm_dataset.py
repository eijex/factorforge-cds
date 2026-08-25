"""Prepare a public-safe candidate-asset manifest for FactorForge-LM research.

This script intentionally does not create a train/validation/test corpus. It only
fingerprints packaged public reference assets so a later release gate can decide
whether a real training-data pipeline is permitted. A local SQLite registration
records this metadata checkpoint; it is not a PostgreSQL or training claim.
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


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_manifest_hash(manifest: dict[str, object]) -> str:
    payload = dict(manifest)
    payload.pop("manifest_hash", None)
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def build_manifest(
    candidate_assets: list[str], data_dir: Path = DATA_DIR
) -> dict[str, object]:
    """Build a byte-pinned manifest without asserting an unperformed split."""
    root = data_dir.resolve()
    fingerprints: list[dict[str, object]] = []
    for relative_path in candidate_assets:
        asset = (root / relative_path).resolve()
        if not asset.is_relative_to(root):
            raise ValueError(f"Candidate asset escapes data directory: {relative_path}")
        if not asset.is_file():
            raise FileNotFoundError(f"Candidate asset does not exist: {asset}")
        fingerprints.append(
            {
                "path": relative_path,
                "size_bytes": asset.stat().st_size,
                "sha256": _sha256_file(asset),
            }
        )

    manifest: dict[str, object] = {
        "manifest_version": "lm-candidate-assets-v1",
        "snapshot_name": "FactorForge-Packaged-Public-Reference-Assets",
        "version": "v1.0",
        "status": "candidate-assets-only-no-training-split",
        "host_scope": "mixed-packaged-reference-assets",
        "sequence_count": None,
        "deduplication_method": None,
        "split_ratio": None,
        "claim_boundary": (
            "This manifest fingerprints packaged candidate assets only. It does not "
            "establish a training corpus, CD-HIT result, train/validation/test split, "
            "held-out benchmark, trained model, or biological-performance claim."
        ),
        "candidate_public_reference_assets": candidate_assets,
        "asset_fingerprints": fingerprints,
    }
    manifest["manifest_hash"] = _canonical_manifest_hash(manifest)
    return manifest


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a byte-pinned candidate-asset manifest for the LM research track."
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT, help="Output JSON path")
    parser.add_argument(
        "--db-path",
        type=Path,
        default=None,
        help="Local SQLite checkpoint path (default: data/db/factorforge_relational.db)",
    )
    parser.add_argument(
        "--skip-db-registration",
        action="store_true",
        help="Write the manifest without registering it in local SQLite",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = build_manifest(collect_candidate_assets())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if not args.skip_db_registration:
        from factorforge.db.connector import FactorForgeDBConnector
        db = FactorForgeDBConnector(db_path=str(args.db_path) if args.db_path else None)
        snapshot_id = db.register_dataset_snapshot(manifest)
        print(
            "Candidate-asset manifest registered in local SQLite "
            f"(snapshot_id={snapshot_id})"
        )

    print(f"Candidate-asset manifest written: {args.out}")
    print(f"Manifest Hash: {manifest['manifest_hash']}")
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
