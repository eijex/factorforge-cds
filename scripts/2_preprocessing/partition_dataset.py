#!/usr/bin/env python3
"""
Split FASTA sequences into train/val/test with stratified TPM bins.
"""

from __future__ import annotations

import argparse
import csv
import logging
import random
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def parse_ratio(ratio: str) -> Tuple[int, int, int]:
    parts = ratio.split(":")
    if len(parts) != 3:
        raise ValueError("Ratio must be in format train:val:test (e.g., 70:10:20)")
    try:
        values = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise ValueError("Ratio values must be integers") from exc
    if sum(values) <= 0:
        raise ValueError("Ratio values must sum to > 0")
    return values  # type: ignore[return-value]


def read_expression(path: Path) -> Dict[str, float]:
    expression: Dict[str, float] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Expression file missing header")
        if "seq_id" not in reader.fieldnames or "TPM" not in reader.fieldnames:
            raise ValueError("Expression file must have columns: seq_id, TPM")
        for row in reader:
            try:
                expression[row["seq_id"]] = float(row["TPM"])
            except ValueError:
                continue
    return expression


def chunk_by_quintile(items: List[Tuple[SeqRecord, float]]) -> List[List[SeqRecord]]:
    if not items:
        return [[] for _ in range(5)]
    items.sort(key=lambda x: x[1])
    total = len(items)
    base = total // 5
    remainder = total % 5
    bins: List[List[SeqRecord]] = []
    idx = 0
    for i in range(5):
        size = base + (1 if i < remainder else 0)
        bin_records = [record for record, _ in items[idx : idx + size]]
        bins.append(bin_records)
        idx += size
    return bins


def split_counts(total: int, ratio: Tuple[int, int, int]) -> Tuple[int, int, int]:
    ratio_sum = sum(ratio)
    train = int(total * ratio[0] / ratio_sum)
    val = int(total * ratio[1] / ratio_sum)
    test = total - train - val
    return train, val, test


def stratified_split(
    records: List[SeqRecord],
    expression: Dict[str, float],
    ratio: Tuple[int, int, int],
    seed: int,
) -> Tuple[List[SeqRecord], List[SeqRecord], List[SeqRecord]]:
    items = [(record, expression.get(record.id, 0.0)) for record in records]
    bins = chunk_by_quintile(items)
    rng = random.Random(seed)

    train_records: List[SeqRecord] = []
    val_records: List[SeqRecord] = []
    test_records: List[SeqRecord] = []

    for bin_records in bins:
        rng.shuffle(bin_records)
        n_train, n_val, n_test = split_counts(len(bin_records), ratio)
        train_records.extend(bin_records[:n_train])
        val_records.extend(bin_records[n_train : n_train + n_val])
        test_records.extend(bin_records[n_train + n_val : n_train + n_val + n_test])

    rng.shuffle(train_records)
    rng.shuffle(val_records)
    rng.shuffle(test_records)
    return train_records, val_records, test_records


def write_fasta(records: List[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "fasta")


def summarize(records: List[SeqRecord]) -> Tuple[int, int, float]:
    total_bp = sum(len(record.seq) for record in records)
    avg_len = total_bp / len(records) if records else 0.0
    return len(records), total_bp, avg_len


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stratified train/val/test split by TPM.")
    parser.add_argument("--input", required=True, help="Input FASTA path.")
    parser.add_argument("--expression", required=True, help="Expression TSV (seq_id, TPM).")
    parser.add_argument("--train", required=True, help="Train FASTA output path.")
    parser.add_argument("--val", required=True, help="Validation FASTA output path.")
    parser.add_argument("--test", required=True, help="Test FASTA output path.")
    parser.add_argument("--ratio", default="70:10:20", help="Split ratio train:val:test.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    try:
        ratio = parse_ratio(args.ratio)
    except ValueError as exc:
        logging.error(str(exc))
        return 1

    try:
        expression = read_expression(Path(args.expression))
    except Exception as exc:
        logging.error("Expression file error: %s", exc)
        return 1

    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        logging.error("No sequences found in %s", args.input)
        return 1

    train_records, val_records, test_records = stratified_split(
        records, expression, ratio, args.seed
    )

    write_fasta(train_records, Path(args.train))
    write_fasta(val_records, Path(args.val))
    write_fasta(test_records, Path(args.test))

    train_stats = summarize(train_records)
    val_stats = summarize(val_records)
    test_stats = summarize(test_records)

    logging.info("Train: %d sequences, %d bp, avg %.1f bp", *train_stats)
    logging.info("Val: %d sequences, %d bp, avg %.1f bp", *val_stats)
    logging.info("Test: %d sequences, %d bp, avg %.1f bp", *test_stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
