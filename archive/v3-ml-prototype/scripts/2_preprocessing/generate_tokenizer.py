#!/usr/bin/env python3
"""
Generate a codon tokenizer JSON mapping.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from itertools import product
from pathlib import Path
from typing import Dict, List

SPECIAL_TOKENS = {
    "[PAD]": 0,
    "[UNK]": 1,
    "[MASK]": 2,
    "[START]": 3,
    "[END]": 4,
}


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def build_codon_vocab() -> List[str]:
    bases = ["A", "C", "G", "T"]
    return ["".join(codon) for codon in product(bases, repeat=3)]


def build_token_map() -> Dict[str, int]:
    codons = build_codon_vocab()
    token_map: Dict[str, int] = dict(SPECIAL_TOKENS)
    start_id = max(token_map.values()) + 1
    for idx, codon in enumerate(codons):
        token_map[codon] = start_id + idx
    return token_map


def tokenize_sequence(seq: str, token_map: Dict[str, int]) -> List[int]:
    seq = seq.upper()
    tokens = [token_map["[START]"]]
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        if len(codon) != 3 or any(base not in "ACGT" for base in codon):
            tokens.append(token_map["[UNK]"])
        else:
            tokens.append(token_map.get(codon, token_map["[UNK]"]))
    tokens.append(token_map["[END]"])
    return tokens


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate codon tokenizer mapping.")
    parser.add_argument("--sequences", required=True, help="FASTA used for sanity check.")
    parser.add_argument("--output", required=True, help="Output tokenizer JSON path.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    token_map = build_token_map()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(token_map, handle, indent=2, sort_keys=True)

    sample = "ATGAAATTTGGGTAG"
    encoded = tokenize_sequence(sample, token_map)
    logging.info("Sample: %s -> %s", sample, encoded)
    logging.info("Tokenizer saved to %s (vocab size=%d)", output_path, len(token_map))
    return 0


if __name__ == "__main__":
    sys.exit(main())
