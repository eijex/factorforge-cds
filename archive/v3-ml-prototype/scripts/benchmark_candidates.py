"""Compare native, v2, v3 fallback, and supplied candidate sequences."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from factorforge.engines import register_builtin_engines
from factorforge.engines.v2.optimizer import RuleBasedOptimizer
from factorforge.engines.v2.utils import parse_fasta_records
from factorforge.engines.v3.pipeline import V3Pipeline
from factorforge.ml.metrics import (
    amino_acid_identity,
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    codon_usage_profile,
    count_internal_stops,
    detect_invalid_codons,
    translate_dna,
)
from factorforge.utils.validation import validate_candidate_sequence


def _load_codon_weights() -> dict[str, float]:
    from factorforge.engines.v3.metrics import load_codon_usage_table

    return load_codon_usage_table().codon_weights


def _read_fasta(path: str | None) -> list[tuple[str, str]]:
    if not path:
        return []
    return parse_fasta_records(Path(path).read_text(encoding="utf-8"))


def _is_dna(sequence: str) -> bool:
    seq = "".join(sequence.upper().split())
    return bool(seq) and len(seq) % 3 == 0 and set(seq) <= set("ATGC")


def _candidate_row(
    name: str,
    dna_sequence: str,
    protein_sequence: str,
    codon_weights: dict[str, float],
) -> dict[str, Any]:
    windows = calculate_gc_windows(dna_sequence)
    window_values = [float(window["gc"]) for window in windows]
    usage = codon_usage_profile(dna_sequence)
    top_codons = ",".join(
        f"{codon}:{int(data['count'])}"
        for codon, data in sorted(
            usage.items(),
            key=lambda item: int(item[1]["count"]),
            reverse=True,
        )[:5]
    )
    validation = validate_candidate_sequence(protein_sequence, dna_sequence)
    return {
        "candidate": name,
        "bp": len(dna_sequence),
        "aa_identity": amino_acid_identity(protein_sequence, dna_sequence),
        "cai": calculate_cai(dna_sequence, codon_weights),
        "gc": calculate_gc(dna_sequence),
        "local_gc_min": min(window_values) if window_values else 0.0,
        "local_gc_max": max(window_values) if window_values else 0.0,
        "internal_stops": count_internal_stops(dna_sequence),
        "invalid_codons": len(detect_invalid_codons(dna_sequence)),
        "validator_passed": validation["passed"],
        "top_codons": top_codons,
    }


def _print_rows(rows: list[dict[str, Any]]) -> None:
    headers = [
        "candidate",
        "bp",
        "CAI",
        "GC%",
        "local GC",
        "AA identity",
        "stops",
        "invalid",
        "valid",
    ]
    rendered = [headers]
    for row in rows:
        rendered.append(
            [
                str(row["candidate"]),
                str(row["bp"]),
                f"{float(row['cai']):.3f}",
                f"{float(row['gc']):.1f}",
                f"{float(row['local_gc_min']):.1f}-{float(row['local_gc_max']):.1f}",
                f"{float(row['aa_identity']) * 100:.1f}%",
                str(row["internal_stops"]),
                str(row["invalid_codons"]),
                str(row["validator_passed"]),
            ]
        )
    widths = [max(len(row[index]) for row in rendered) for index in range(len(headers))]
    for row_index, row in enumerate(rendered):
        print(" | ".join(value.ljust(widths[index]) for index, value in enumerate(row)))
        if row_index == 0:
            print("-|-".join("-" * width for width in widths))


def _write_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark FactorForge candidate CDS sequences.")
    parser.add_argument("--protein", help="Protein sequence.")
    parser.add_argument("--protein-fasta", help="Protein FASTA; first record is used.")
    parser.add_argument("--native-fasta", help="Native CDS FASTA; all records are compared.")
    parser.add_argument("--candidate-fasta", help="Candidate CDS FASTA; all records are compared.")
    parser.add_argument("--include-v3", action="store_true", help="Include v3 no-checkpoint fallback.")
    parser.add_argument("--csv", nargs="?", const="experiments/results/baselines/candidate_metrics.csv")
    args = parser.parse_args()

    if args.protein:
        protein_sequence = "".join(args.protein.upper().split())
    elif args.protein_fasta:
        protein_sequence = _read_fasta(args.protein_fasta)[0][1]
    elif args.native_fasta:
        native = _read_fasta(args.native_fasta)[0][1]
        protein_sequence = translate_dna(native).rstrip("*")
    else:
        raise ValueError("Provide --protein, --protein-fasta, or --native-fasta")

    codon_weights = _load_codon_weights()
    candidates: list[tuple[str, str]] = []
    for name, sequence in _read_fasta(args.native_fasta):
        if _is_dna(sequence):
            candidates.append((f"native:{name}", sequence))

    try:
        v2_result = RuleBasedOptimizer().optimize(protein_sequence, profile="balanced")
        candidates.append(("v2", v2_result.sequence))
    except Exception as exc:
        print(f"v2 unavailable: {exc}", file=sys.stderr)

    if args.include_v3:
        try:
            register_builtin_engines()
            v3_result = V3Pipeline().run(protein_sequence)
            candidates.append(("v3_fallback", v3_result.sequence))
        except Exception as exc:
            print(f"v3 unavailable: {exc}", file=sys.stderr)

    for name, sequence in _read_fasta(args.candidate_fasta):
        if _is_dna(sequence):
            candidates.append((f"candidate:{name}", sequence))

    rows = [_candidate_row(name, sequence, protein_sequence, codon_weights) for name, sequence in candidates]
    if not rows:
        raise ValueError("No candidate DNA sequences available for benchmarking")

    _print_rows(rows)
    if args.csv:
        output_path = ROOT / args.csv if not Path(args.csv).is_absolute() else Path(args.csv)
        _write_csv(rows, output_path)
        print(f"\nCSV written: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

