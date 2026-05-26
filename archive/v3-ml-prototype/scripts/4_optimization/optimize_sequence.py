#!/usr/bin/env python3
"""
Optimize a protein sequence using a trained ML codon optimizer.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple

from Bio import SeqIO

from factorforge.engines.ml.plant_optimizer import PlantCodonOptimizer


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def read_protein_input(value: str) -> Tuple[str, str]:
    path = Path(value)
    if path.exists():
        records = list(SeqIO.parse(str(path), "fasta"))
        if not records:
            raise ValueError(f"No FASTA records found in {path}")
        record = records[0]
        return record.id, str(record.seq)
    return "Protein", value


def write_fasta(sequence: str, path: Path, name: str, description: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{name} {description}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i : i + 60] + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optimize a protein sequence.")
    parser.add_argument("--input", required=True, help="Protein FASTA file or sequence string.")
    parser.add_argument("--model", required=True, help="Path to trained model checkpoint.")
    parser.add_argument("--tokenizer", required=True, help="Path to tokenizer.json.")
    parser.add_argument("--codon-table", required=True, help="Path to codon frequency CSV.")
    parser.add_argument("--organism", default="nbenthamiana", help="Organism name.")
    parser.add_argument("--output", required=True, help="Output FASTA file.")
    parser.add_argument("--report", default=None, help="Optional report file path.")
    parser.add_argument("--beam-size", type=int, default=5, help="Beam search size.")
    parser.add_argument("--original-dna", default=None, help="Optional original DNA sequence.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    try:
        protein_name, protein_seq = read_protein_input(args.input)
    except Exception as exc:
        logging.error("Input error: %s", exc)
        return 1

    logging.info("Loaded protein: %s (%d aa)", protein_name, len(protein_seq))
    optimizer = PlantCodonOptimizer(
        model_path=args.model,
        codon_table_path=args.codon_table,
        tokenizer_path=args.tokenizer,
        organism=args.organism,
    )

    try:
        logging.info("Optimizing with beam size %d", args.beam_size)
        optimized_dna = optimizer.optimize(protein_seq, beam_size=args.beam_size)
    except Exception as exc:
        logging.error("Optimization failed: %s", exc)
        return 1

    output_path = Path(args.output)
    description = f"optimized_for={args.organism}"
    write_fasta(optimized_dna, output_path, protein_name, description)
    logging.info("Optimized FASTA written to %s", output_path)

    original_dna = args.original_dna
    if original_dna:
        original_dna = original_dna.strip().upper()
        if len(original_dna) % 3 != 0:
            logging.error("Original DNA length must be divisible by 3.")
            return 1

    if args.report:
        report_path = Path(args.report)
        if original_dna is None:
            original_dna = optimizer._reverse_translate(protein_seq)
        optimizer.generate_report(
            protein_name=protein_name,
            protein_seq=protein_seq,
            original_dna=original_dna,
            optimized_dna=optimized_dna,
            output_path=str(report_path),
        )
        logging.info("Report written to %s", report_path)
    elif original_dna:
        metrics = optimizer.compare_sequences(original_dna, optimized_dna)
        logging.info(
            "Change rate: %.2f%% | CAI: %.4f -> %.4f",
            metrics["change_rate"],
            metrics["original_cai"],
            metrics["optimized_cai"],
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
