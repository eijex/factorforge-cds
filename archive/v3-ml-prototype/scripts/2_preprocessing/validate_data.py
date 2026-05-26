#!/usr/bin/env python3
"""
Validate train/val/test FASTA splits.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

from Bio import SeqIO

STOP_CODONS = {"TAA", "TAG", "TGA"}


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def load_sequences(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    for record in SeqIO.parse(path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


def gc_content(seqs: List[str]) -> float:
    total = sum(len(seq) for seq in seqs)
    if total == 0:
        return 0.0
    gc = sum(seq.count("G") + seq.count("C") for seq in seqs)
    return (gc / total) * 100.0


def count_stops(seqs: List[str]) -> Tuple[Dict[str, int], int]:
    distribution = {codon: 0 for codon in STOP_CODONS}
    premature = 0
    for seq in seqs:
        for i in range(0, len(seq) - 3, 3):
            codon = seq[i : i + 3]
            if codon in STOP_CODONS:
                distribution[codon] += 1
                premature += 1
    return distribution, premature


def validate_lengths(seqs: List[str]) -> int:
    return sum(1 for seq in seqs if len(seq) % 3 != 0)


def report_split(name: str, sequences: Dict[str, str]) -> None:
    seq_list = list(sequences.values())
    total_bp = sum(len(seq) for seq in seq_list)
    total_codons = total_bp // 3
    avg_len = total_bp / len(seq_list) if seq_list else 0.0
    gc = gc_content(seq_list)
    stop_dist, premature = count_stops(seq_list)
    invalid_len = validate_lengths(seq_list)

    print(f"{name}: {len(seq_list)} sequences, {total_codons} codons, GC={gc:.1f}%")
    print(f"  Avg length: {avg_len:.1f} bp, Invalid length: {invalid_len}")
    print(
        f"  Premature stops: {premature} (TAA={stop_dist['TAA']}, "
        f"TAG={stop_dist['TAG']}, TGA={stop_dist['TGA']})"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate train/val/test datasets.")
    parser.add_argument("--train", required=True, help="Train FASTA path.")
    parser.add_argument("--val", required=True, help="Validation FASTA path.")
    parser.add_argument("--test", required=True, help="Test FASTA path.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    train = load_sequences(Path(args.train))
    val = load_sequences(Path(args.val))
    test = load_sequences(Path(args.test))

    print("=== Dataset Validation ===")
    report_split("Train", train)
    report_split("Val", val)
    report_split("Test", test)

    ids_train = set(train.keys())
    ids_val = set(val.keys())
    ids_test = set(test.keys())

    overlap_ids = (ids_train & ids_val) | (ids_train & ids_test) | (ids_val & ids_test)
    seqs_train = set(train.values())
    seqs_val = set(val.values())
    seqs_test = set(test.values())
    overlap_seqs = (seqs_train & seqs_val) | (seqs_train & seqs_test) | (seqs_val & seqs_test)

    issues = 0
    if overlap_ids:
        print(f"Overlap IDs detected: {len(overlap_ids)}")
        issues += 1
    if overlap_seqs:
        print(f"Overlap sequences detected: {len(overlap_seqs)}")
        issues += 1

    invalid_lengths = (
        validate_lengths(list(train.values()))
        + validate_lengths(list(val.values()))
        + validate_lengths(list(test.values()))
    )
    if invalid_lengths:
        print(f"Invalid length sequences: {invalid_lengths}")
        issues += 1

    _, premature_train = count_stops(list(train.values()))
    _, premature_val = count_stops(list(val.values()))
    _, premature_test = count_stops(list(test.values()))
    premature_total = premature_train + premature_val + premature_test
    if premature_total:
        print(f"Premature stop codons detected: {premature_total}")
        issues += 1

    if issues == 0:
        print("No overlaps detected \u2713")
        print("All sequences valid \u2713")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
