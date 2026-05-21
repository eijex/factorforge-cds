"""Diagnose v2/v3 sequence metrics without training.

This script accepts a protein sequence, DNA sequence, or FASTA file and prints a
small comparison table for native, v2, and v3 fallback candidates when available.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from factorforge.engines import register_builtin_engines
from factorforge.engines.v2.optimizer import RuleBasedOptimizer
from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile, ReverseTranslator
from factorforge.engines.v2.utils import parse_fasta_records
from factorforge.engines.v3.metrics import compute_cai, compute_gc, load_codon_usage_table
from factorforge.engines.v3.pipeline import V3Pipeline


STOP_CODONS = {"TAA", "TAG", "TGA"}


def _load_sequence(args: argparse.Namespace) -> tuple[str, str]:
    if args.protein:
        return "input", args.protein.strip().upper()
    if args.sequence:
        return "input", args.sequence.strip().upper()
    if args.fasta:
        content = Path(args.fasta).read_text(encoding="utf-8")
        records = parse_fasta_records(content)
        if not records:
            raise ValueError("No FASTA records found")
        return records[0]
    raise ValueError("Provide --protein, --sequence, or --fasta")


def _is_dna(sequence: str) -> bool:
    return bool(sequence) and set(sequence.upper()) <= set("ATGC") and len(sequence) % 3 == 0


def _codons(sequence: str) -> list[str]:
    seq = sequence.upper().replace("U", "T")
    return [seq[index : index + 3] for index in range(0, len(seq) - len(seq) % 3, 3)]


def _translate(dna_sequence: str, codon_to_aa: dict[str, str]) -> str:
    return "".join(codon_to_aa.get(codon, "X") for codon in _codons(dna_sequence))


def _local_gc_windows(sequence: str, window_size: int = 60, step: int = 30) -> list[float]:
    seq = sequence.upper()
    if len(seq) < window_size:
        return [compute_gc(seq)] if seq else []
    return [
        compute_gc(seq[start : start + window_size])
        for start in range(0, len(seq) - window_size + 1, step)
    ]


def _summarize_candidate(
    name: str,
    dna_sequence: str,
    protein_sequence: str,
    codon_to_aa: dict[str, str],
    codon_weights_source,
) -> dict[str, object]:
    translated = _translate(dna_sequence, codon_to_aa)
    comparable_length = min(len(protein_sequence), len(translated))
    matches = sum(
        1
        for expected, observed in zip(
            protein_sequence[:comparable_length],
            translated[:comparable_length],
        )
        if expected == observed
    )
    identity = matches / len(protein_sequence) if protein_sequence else 0.0
    internal_codons = _codons(dna_sequence)[:-1]
    invalid_codons = [codon for codon in _codons(dna_sequence) if codon not in codon_to_aa]
    stop_count = sum(1 for codon in internal_codons if codon in STOP_CODONS)
    windows = _local_gc_windows(dna_sequence)
    usage = Counter(_codons(dna_sequence))
    most_common = ",".join(f"{codon}:{count}" for codon, count in usage.most_common(5))
    return {
        "name": name,
        "length": len(dna_sequence),
        "cai": compute_cai(dna_sequence, codon_weights_source),
        "gc": compute_gc(dna_sequence),
        "local_gc_min": min(windows) if windows else 0.0,
        "local_gc_max": max(windows) if windows else 0.0,
        "aa_identity": identity,
        "internal_stops": stop_count,
        "invalid_codons": len(invalid_codons),
        "top_codons": most_common,
    }


def _print_table(rows: Iterable[dict[str, object]]) -> None:
    headers = [
        "candidate",
        "bp",
        "CAI",
        "GC%",
        "local GC",
        "AA identity",
        "stops",
        "invalid",
        "top codons",
    ]
    rendered = [headers]
    for row in rows:
        rendered.append(
            [
                str(row["name"]),
                str(row["length"]),
                f"{float(row['cai']):.3f}",
                f"{float(row['gc']):.1f}",
                f"{float(row['local_gc_min']):.1f}-{float(row['local_gc_max']):.1f}",
                f"{float(row['aa_identity']) * 100:.1f}%",
                str(row["internal_stops"]),
                str(row["invalid_codons"]),
                str(row["top_codons"]),
            ]
        )

    widths = [max(len(row[index]) for row in rendered) for index in range(len(headers))]
    for row_index, row in enumerate(rendered):
        print(" | ".join(value.ljust(widths[index]) for index, value in enumerate(row)))
        if row_index == 0:
            print("-|-".join("-" * width for width in widths))


def _run_v2(protein_sequence: str, profile: str, force_protein: bool) -> str:
    if force_protein:
        translator = ReverseTranslator()
        opt_profile = OptimizationProfile(profile)
        return translator.generate_candidates(protein_sequence, profile=opt_profile, n=1)[0]["sequence"]

    result = RuleBasedOptimizer().optimize(protein_sequence, profile=profile)
    return result.sequence


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose FactorForge v2/v3 metrics.")
    parser.add_argument("--protein", help="Protein amino acid sequence.")
    parser.add_argument("--sequence", help="Protein or DNA sequence.")
    parser.add_argument("--fasta", help="Path to FASTA file; first record is used.")
    parser.add_argument("--profile", default="balanced", help="v2 optimization profile.")
    args = parser.parse_args()

    record_id, input_sequence = _load_sequence(args)
    table = load_codon_usage_table()
    protein_sequence = (
        _translate(input_sequence, table.codon_to_aa).rstrip("*")
        if _is_dna(input_sequence)
        else input_sequence
    )

    rows: list[dict[str, object]] = []
    if _is_dna(input_sequence):
        rows.append(
            _summarize_candidate("native", input_sequence, protein_sequence, table.codon_to_aa, table)
        )

    try:
        v2_sequence = _run_v2(protein_sequence, args.profile, force_protein=bool(args.protein))
        rows.append(_summarize_candidate("v2", v2_sequence, protein_sequence, table.codon_to_aa, table))
    except Exception as exc:
        print(f"v2 unavailable: {exc}", file=sys.stderr)

    try:
        register_builtin_engines()
        v3_result = V3Pipeline().run(protein_sequence)
        rows.append(_summarize_candidate("v3_fallback", v3_result.sequence, protein_sequence, table.codon_to_aa, table))
    except Exception as exc:
        print(f"v3 unavailable: {exc}", file=sys.stderr)

    print(f"Record: {record_id}")
    print(f"Protein length: {len(protein_sequence)} aa")
    _print_table(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
