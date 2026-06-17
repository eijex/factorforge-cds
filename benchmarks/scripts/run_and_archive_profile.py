"""Run the formal benchmark for one codon source profile and archive its summary."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUMMARY_PATH = ROOT / "benchmarks" / "results" / "v3.2.0" / "benchmark_summary.json"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-id", required=True)
    parser.add_argument("--codon-table-path", default=None)
    parser.add_argument("--manifest-path", default=None)
    parser.add_argument("--seed", type=int, default=320)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "reproducibility" / "benchmark_v3.2.2")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    command = [
        sys.executable,
        str(ROOT / "benchmarks" / "run_benchmark.py"),
        "--dataset",
        "formal",
        "--mode",
        "formal",
        "--seed",
        str(args.seed),
        "--source-profile-id",
        args.profile_id,
    ]
    if args.codon_table_path:
        command.extend(["--codon-table-path", args.codon_table_path])
        if args.manifest_path:
            command.extend(["--source-profile-manifest", args.manifest_path])

    completed = subprocess.run(command, cwd=ROOT)
    if completed.returncode != 0:
        return completed.returncode
    if not SUMMARY_PATH.exists():
        raise SystemExit(f"Missing benchmark summary after run: {SUMMARY_PATH}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    destination = args.out_dir / f"benchmark_summary_profile_{args.profile_id}.frozen.json"
    shutil.copy2(SUMMARY_PATH, destination)
    print(f"Archived: {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
