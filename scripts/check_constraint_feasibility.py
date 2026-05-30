"""
check_constraint_feasibility.py — DP feasibility diagnostic tool

Computes the achievable CAI/GC% range for a protein sequence using the
DP Feasibility Engine, without running full optimization.

Usage:
    python scripts/check_constraint_feasibility.py "MVSK..."
    python scripts/check_constraint_feasibility.py --gc-min 55 --gc-max 65 "MVSK..."
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from factorforge.analysis.feasibility import analyze_feasibility  # noqa: E402
from factorforge.analysis.metrics import load_codon_usage_table  # noqa: E402
from factorforge.engines.profile.utils import parse_fasta_records  # noqa: E402


def _load_protein(args: argparse.Namespace) -> str:
    if args.protein:
        return "".join(args.protein.upper().split())
    if args.fasta:
        records = parse_fasta_records(Path(args.fasta).read_text(encoding="utf-8"))
        return records[0][1]
    raise ValueError("Provide --protein or --fasta")


def main() -> int:
    parser = argparse.ArgumentParser(description="Check synonymous CAI/GC feasibility.")
    parser.add_argument("--protein", help="Protein sequence.")
    parser.add_argument("--fasta", help="Protein FASTA path; first record is used.")
    parser.add_argument("--codon-table", help="Optional codon table JSON path.")
    parser.add_argument("--cai-target", type=float, default=0.92)
    parser.add_argument("--gc-low", type=float, default=0.41)
    parser.add_argument("--gc-high", type=float, default=0.44)
    args = parser.parse_args()

    table = load_codon_usage_table(Path(args.codon_table) if args.codon_table else None)
    protein = _load_protein(args)
    result = analyze_feasibility(
        protein,
        table.codon_weights,
        target_cai=args.cai_target,
        target_gc_low=args.gc_low,
        target_gc_high=args.gc_high,
    )

    target = result["target"]
    print(f"Protein length: {result['protein_length']} aa")
    print(
        "Possible GC range: "
        f"{result['minimum_possible_gc']:.1f}% - {result['maximum_possible_gc']:.1f}%"
    )
    print(
        "Max CAI without GC constraint: "
        f"{result['maximum_achievable_cai_without_gc']:.3f}"
    )
    print(
        "Target feasible: "
        f"{target['feasible']} "
        f"(CAI >= {target['cai']:.3f}, GC {target['gc_low']:.1f}-{target['gc_high']:.1f}%)"
    )
    print("\nGC range   | feasible | max CAI | best GC")
    print("-----------|----------|---------|--------")
    for label, data in result["ranges"].items():
        candidate = data["best_candidate"]
        max_cai = "-" if data["max_cai"] is None else f"{data['max_cai']:.3f}"
        best_gc = "-" if candidate is None else f"{candidate['gc']:.1f}%"
        print(f"{label.ljust(10)} | {str(data['feasible']).ljust(8)} | {max_cai.ljust(7)} | {best_gc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
