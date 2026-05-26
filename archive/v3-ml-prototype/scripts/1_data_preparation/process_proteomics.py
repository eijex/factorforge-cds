#!/usr/bin/env python3
"""
Process PRIDE proteomics data (PXD042916) and export filtered FASTA.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, Iterable, List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def guess_delimiter(path: Path, delimiter: str | None) -> str:
    if delimiter:
        return delimiter
    if path.suffix.lower() == ".csv":
        return ","
    return "\t"


def parse_rows(
    path: Path,
    delimiter: str,
    psm_col: str,
    fdr_col: str,
    seq_col: str,
    id_col: str | None,
    min_psm: float,
    max_fdr: float,
) -> List[SeqRecord]:
    records: List[SeqRecord] = []
    seen = set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for idx, row in enumerate(reader, start=1):
            try:
                psm = float(row.get(psm_col, "0") or 0)
                fdr = float(row.get(fdr_col, "1") or 1)
            except ValueError:
                continue
            if psm <= min_psm or fdr >= max_fdr:
                continue
            sequence = (row.get(seq_col) or "").strip()
            if not sequence:
                continue
            sequence = sequence.replace(" ", "").replace("\t", "").upper()
            if sequence in seen:
                continue
            seen.add(sequence)
            record_id = (row.get(id_col) or f"protein_{idx}") if id_col else f"protein_{idx}"
            records.append(SeqRecord(Seq(sequence), id=record_id, description=""))
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter PRIDE proteomics data and export FASTA.")
    parser.add_argument("--input", required=True, help="Input PSM table (TSV/CSV).")
    parser.add_argument("--output", required=True, help="Output FASTA path.")
    parser.add_argument("--psm-column", default="PSM", help="Column name for PSM count.")
    parser.add_argument("--fdr-column", default="FDR", help="Column name for FDR.")
    parser.add_argument("--sequence-column", default="Sequence", help="Column name for peptide.")
    parser.add_argument("--id-column", default="Protein", help="Column name for protein ID.")
    parser.add_argument("--min-psm", type=float, default=2.0, help="Minimum PSM threshold.")
    parser.add_argument("--max-fdr", type=float, default=0.01, help="Maximum FDR threshold.")
    parser.add_argument("--delimiter", default=None, help="Override delimiter.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    input_path = Path(args.input)
    output_path = Path(args.output)
    delimiter = guess_delimiter(input_path, args.delimiter)

    records = parse_rows(
        input_path,
        delimiter,
        args.psm_column,
        args.fdr_column,
        args.sequence_column,
        args.id_column,
        args.min_psm,
        args.max_fdr,
    )

    if not records:
        logging.error("No records passed filters. Check column names and thresholds.")
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output_path, "fasta")
    logging.info("Wrote %d sequences to %s", len(records), output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
