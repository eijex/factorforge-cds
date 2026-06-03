"""Run a reproducible FactorForge benchmark."""

from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from factorforge.engines.profile.optimizer import RuleBasedOptimizer  # noqa: E402
from factorforge.utils.sequence_validator import validate_cds_output  # noqa: E402


DEFAULT_BENCHMARK_SIZE = 3876
DEFAULT_SEQUENCE_LENGTH = 100
PANEL_SEED = 123
OPTIMIZER_SEED = 123
GC_ENRICHMENT_RATE = 0.08
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
GC_BALANCE_AA = "APGWR"
VALID_PROFILES = ("balanced", "high_cai", "gc_target")


@dataclass(frozen=True)
class BenchmarkRow:
    sequence_id: str
    profile: str
    cai: float
    gc_percent: float
    aa_identity: float
    validator_passed: bool
    output_length_nt: int


def build_panel(size: int, length: int = DEFAULT_SEQUENCE_LENGTH) -> list[tuple[str, str]]:
    """Build a deterministic plant protein panel from the codon reference alphabet."""
    rng = random.Random(PANEL_SEED)
    panel: list[tuple[str, str]] = []

    for index in range(size):
        residues = ["M"]
        for _ in range(length - 1):
            alphabet = GC_BALANCE_AA if rng.random() < GC_ENRICHMENT_RATE else STANDARD_AA
            residues.append(rng.choice(alphabet))
        panel.append((f"benchmark_{index + 1:04d}", "".join(residues)))

    return panel


def run_profile(profile: str, panel: list[tuple[str, str]]) -> list[BenchmarkRow]:
    random.seed(OPTIMIZER_SEED)
    optimizer = RuleBasedOptimizer()
    rows: list[BenchmarkRow] = []

    for sequence_id, protein in panel:
        result = optimizer.optimize(protein, profile=profile, scan_mode="fast")
        validation = validate_cds_output(protein, result.sequence)
        rows.append(
            BenchmarkRow(
                sequence_id=sequence_id,
                profile=profile,
                cai=float(result.metrics["cai"]),
                gc_percent=float(result.metrics.get("gc_percent", result.metrics["gc_content"])),
                aa_identity=float(validation["aa_identity"]),
                validator_passed=bool(validation["passed"]),
                output_length_nt=len(result.sequence),
            )
        )

    return rows


def summarize(rows: list[BenchmarkRow]) -> dict[str, float]:
    if not rows:
        raise ValueError("Benchmark requires at least one sequence")

    return {
        "sequences": float(len(rows)),
        "cai_mean": mean(row.cai for row in rows),
        "gc_mean": mean(row.gc_percent for row in rows),
        "aa_identity": mean(row.aa_identity for row in rows) * 100.0,
        "validator_pass_rate": mean(1.0 if row.validator_passed else 0.0 for row in rows)
        * 100.0,
    }


def print_summary(profile: str, rows: list[BenchmarkRow], runtime: float) -> None:
    summary = summarize(rows)
    print(f"FactorForge Benchmark - {profile} profile (N. benthamiana)")
    print(f"Sequences: {int(summary['sequences'])}")
    print(f"CAI mean:  {summary['cai_mean']:.3f}  (target >= 0.75)")
    print(f"GC% mean:  {summary['gc_mean']:.2f}% (target 55-65%)")
    print(f"AA identity: {summary['aa_identity']:.0f}%")
    print(f"Validator pass rate: {summary['validator_pass_rate']:.0f}%")
    print(f"Runtime: {runtime:.1f}s")


def write_csv(path: Path, rows: list[BenchmarkRow]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_id",
                "profile",
                "cai",
                "gc_percent",
                "aa_identity",
                "validator_passed",
                "output_length_nt",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sequence_id": row.sequence_id,
                    "profile": row.profile,
                    "cai": f"{row.cai:.3f}",
                    "gc_percent": f"{row.gc_percent:.2f}",
                    "aa_identity": f"{row.aa_identity:.4f}",
                    "validator_passed": row.validator_passed,
                    "output_length_nt": row.output_length_nt,
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the FactorForge reproducible benchmark.")
    parser.add_argument("--n", type=int, default=DEFAULT_BENCHMARK_SIZE, help="Number of sequences")
    parser.add_argument("--profile", choices=VALID_PROFILES, default="balanced")
    parser.add_argument("--output", type=Path, help="Optional CSV output path")
    parser.add_argument(
        "--compare-profiles",
        help="Comma-separated profiles to benchmark, e.g. balanced,high_cai,gc_target",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.n < 1:
        raise ValueError("--n must be >= 1")

    profiles = [args.profile]
    if args.compare_profiles:
        profiles = [profile.strip() for profile in args.compare_profiles.split(",") if profile.strip()]
        invalid = sorted(set(profiles) - set(VALID_PROFILES))
        if invalid:
            raise ValueError(f"Unsupported profiles: {', '.join(invalid)}")

    panel = build_panel(args.n)
    all_rows: list[BenchmarkRow] = []

    for index, profile in enumerate(profiles):
        started = time.perf_counter()
        rows = run_profile(profile, panel)
        runtime = time.perf_counter() - started
        if index:
            print()
        print_summary(profile, rows, runtime)
        all_rows.extend(rows)

    if args.output:
        write_csv(args.output, all_rows)
        print(f"\nCSV written: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
