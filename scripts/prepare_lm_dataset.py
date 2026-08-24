"""Prepare a public-safe manifest for the experimental FactorForge-LM track.

This script intentionally does not create a train/validation/test corpus. It only
summarizes candidate public reference assets so a later release gate can decide
whether a real training-data pipeline is permitted and how split/provenance hashes
should be recorded.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "src" / "factorforge" / "data"
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
    """Build a non-training manifest for review and release gating."""
    return {
        "manifest_version": "lm-candidate-assets-v0",
        "status": "candidate-assets-only-no-training-split",
        "claim_boundary": (
            "This manifest does not establish a training corpus, CD-HIT split, "
            "model benchmark, or biological-performance claim."
        ),
        "candidate_public_reference_assets": candidate_assets,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a candidate-asset manifest for the experimental LM track."
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT, help="Output JSON path")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = build_manifest(collect_candidate_assets())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Candidate manifest written: {args.out}")
    print("Status: candidate assets only; no training split or benchmark created")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
