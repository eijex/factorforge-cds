"""Filter SGN CDS records by CAI and write v3 Run 2 training pairs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from factorforge.engines.v3.metrics import compute_cai, load_codon_usage_table


DEFAULT_INPUT = Path("data/raw/sgn_nbenthamiana_cds.fasta")
DEFAULT_OUTPUT = Path("data/training/training_pairs_v3_run2.jsonl")
DEFAULT_CODON_TABLE = Path("data/nbenthamiana_codons.json")
STOP_CODONS = {"TAA", "TAG", "TGA"}
STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]


def strip_terminal_stop(cds: str) -> str:
    seq = cds.upper().replace("U", "T")
    if seq[-3:] in STOP_CODONS:
        return seq[:-3]
    return seq


def translate_cds(cds_without_stop: str) -> str | None:
    amino_acids: list[str] = []
    for index in range(0, len(cds_without_stop), 3):
        codon = cds_without_stop[index : index + 3]
        aa = STANDARD_TABLE.forward_table.get(codon)
        if aa is None:
            return None
        amino_acids.append(aa)
    return "".join(amino_acids)


def filter_pairs(
    input_path: Path,
    output_path: Path | None,
    threshold: float,
    codon_table_path: Path,
) -> int:
    table = load_codon_usage_table(codon_table_path)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    handle = output_path.open("w", encoding="utf-8") if output_path is not None else None
    try:
        for record in SeqIO.parse(input_path, "fasta"):
            row = pair_from_record(record, table, threshold)
            if row is None:
                continue
            if handle is not None:
                handle.write(json.dumps(row, ensure_ascii=False) + "\n")
            written += 1
    finally:
        if handle is not None:
            handle.close()

    return written


def pair_from_record(record, table, threshold: float) -> dict[str, str] | None:
    cds = str(record.seq).upper().replace("U", "T")
    coding_cds = strip_terminal_stop(cds)
    if len(coding_cds) == 0 or len(coding_cds) % 3 != 0:
        return None
    cai = compute_cai(coding_cds, table)
    if cai < threshold:
        return None
    aa_seq = translate_cds(coding_cds)
    if not aa_seq:
        return None
    return {
        "id": record.id,
        "protein_id": record.id,
        "aa_seq": aa_seq,
        "sequence": aa_seq,
        "cds": coding_cds,
        "codon_sequence": " ".join(
            coding_cds[index : index + 3] for index in range(0, len(coding_cds), 3)
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter CDS FASTA by CAI.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--threshold", type=float, default=0.55)
    parser.add_argument("--codon-table", type=Path, default=DEFAULT_CODON_TABLE)
    parser.add_argument("--count", action="store_true", help="Print filtered pair count only.")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Input FASTA not found: {args.input}")
        return 0 if args.count else 1

    output_path = None if args.count else args.output
    count = filter_pairs(args.input, output_path, args.threshold, args.codon_table)
    if args.count:
        print(count)
    else:
        print(f"Saved {count} CAI-filtered pairs to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
