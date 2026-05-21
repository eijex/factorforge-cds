#!/usr/bin/env python3
"""
Build codon frequency and RSCU table from CDS and expression data.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple

from Bio import SeqIO
from Bio.Data import CodonTable


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def read_expression(path: Path, id_col: str, expr_col: str) -> Dict[str, float]:
    expression: Dict[str, float] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if id_col not in reader.fieldnames or expr_col not in reader.fieldnames:
            raise ValueError(f"Missing columns in expression file: {id_col}, {expr_col}")
        for row in reader:
            try:
                expression[row[id_col]] = float(row[expr_col])
            except ValueError:
                continue
    return expression


def calculate_codon_counts(
    cds_path: Path,
    expression: Dict[str, float],
    default_tpm: float,
) -> Dict[str, float]:
    counts: Dict[str, float] = defaultdict(float)
    for record in SeqIO.parse(cds_path, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < 3:
            continue
        if len(seq) % 3 != 0:
            logging.warning("Skipping %s: length not divisible by 3", record.id)
            continue
        weight = expression.get(record.id, default_tpm)
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i : i + 3]
            if "N" in codon or "-" in codon:
                continue
            counts[codon] += weight
    return counts


def rscu_table(counts: Dict[str, float]) -> Dict[str, Tuple[str, float]]:
    table = CodonTable.unambiguous_dna_by_id[1]
    amino_acid_for = {**table.forward_table}
    for stop in table.stop_codons:
        amino_acid_for[stop] = "*"

    aa_codons: Dict[str, list] = defaultdict(list)
    for codon, aa in amino_acid_for.items():
        if aa == "*":
            continue
        aa_codons[aa].append(codon)

    rscu: Dict[str, Tuple[str, float]] = {}
    for aa, codons in aa_codons.items():
        total = sum(counts.get(codon, 0.0) for codon in codons)
        if total == 0:
            for codon in codons:
                rscu[codon] = (aa, 0.0)
            continue
        expected = total / len(codons)
        for codon in codons:
            observed = counts.get(codon, 0.0)
            rscu[codon] = (aa, observed / expected if expected else 0.0)
    return rscu


def write_table(
    counts: Dict[str, float],
    rscu: Dict[str, Tuple[str, float]],
    output_path: Path,
) -> None:
    total_counts = sum(counts.values()) or 1.0
    rows = []
    for codon, (aa, rscu_value) in rscu.items():
        count = counts.get(codon, 0.0)
        rows.append(
            (
                codon,
                aa,
                count,
                count / total_counts,
                rscu_value,
            )
        )

    rows.sort(key=lambda x: (x[1], x[0]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Codon", "AminoAcid", "Count", "Frequency", "RSCU"])
        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build codon frequency and RSCU table.")
    parser.add_argument("--cds", required=True, help="CDS FASTA path.")
    parser.add_argument("--expression", required=True, help="Expression table (quant.sf).")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument("--id-column", default="Name", help="Expression ID column.")
    parser.add_argument("--expr-column", default="TPM", help="Expression value column.")
    parser.add_argument("--default-tpm", type=float, default=1.0, help="Default TPM.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    cds_path = Path(args.cds)
    expr_path = Path(args.expression)
    output_path = Path(args.output)

    try:
        expression = read_expression(expr_path, args.id_column, args.expr_column)
    except Exception as exc:
        logging.error("Failed to read expression data: %s", exc)
        return 1

    counts = calculate_codon_counts(cds_path, expression, args.default_tpm)
    if not counts:
        logging.error("No codon counts produced.")
        return 1

    rscu = rscu_table(counts)
    write_table(counts, rscu, output_path)
    logging.info("Codon table written to %s", output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
