#!/usr/bin/env python3
"""
Run tblastn to match proteins against a genome and extract CDS sequences.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found on PATH: {tool}")


def run_command(cmd: List[str]) -> None:
    logging.debug("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def blast_db_exists(prefix: Path) -> bool:
    return any((prefix.with_suffix(suf)).exists() for suf in [".nhr", ".nin", ".nsq"])


def make_blast_db(genome: Path, db_prefix: Path) -> None:
    if blast_db_exists(db_prefix):
        logging.info("BLAST DB exists: %s", db_prefix)
        return
    run_command(["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db_prefix)])


def run_tblastn(
    proteins: Path,
    db_prefix: Path,
    out_table: Path,
    evalue: float,
    threads: int,
) -> None:
    run_command(
        [
            "tblastn",
            "-query",
            str(proteins),
            "-db",
            str(db_prefix),
            "-evalue",
            str(evalue),
            "-outfmt",
            "6 qseqid sseqid sstart send evalue bitscore",
            "-num_threads",
            str(threads),
            "-out",
            str(out_table),
        ]
    )


def parse_hits(out_table: Path, evalue: float) -> List[Tuple[str, str, int, int, float, float]]:
    hits: List[Tuple[str, str, int, int, float, float]] = []
    with out_table.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            qseqid, sseqid, sstart, send, hit_eval, bitscore = parts
            hit_eval_f = float(hit_eval)
            if hit_eval_f > evalue:
                continue
            hits.append((qseqid, sseqid, int(sstart), int(send), hit_eval_f, float(bitscore)))
    return hits


def extract_cds(
    hits: List[Tuple[str, str, int, int, float, float]],
    genome_records: Dict[str, SeqRecord],
) -> List[SeqRecord]:
    extracted: List[SeqRecord] = []
    for idx, (qseqid, sseqid, sstart, send, hit_eval, bitscore) in enumerate(hits, start=1):
        record = genome_records.get(sseqid)
        if record is None:
            continue
        start = min(sstart, send) - 1
        end = max(sstart, send)
        seq = record.seq[start:end]
        if sstart > send:
            seq = seq.reverse_complement()
        record_id = f"{qseqid}|{sseqid}|hit{idx}"
        description = f"pos={sstart}-{send} evalue={hit_eval:.2e} bitscore={bitscore:.2f}"
        extracted.append(SeqRecord(Seq(str(seq)), id=record_id, description=description))
    return extracted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run tblastn and extract CDS sequences.")
    parser.add_argument("--proteins", required=True, help="Protein FASTA path.")
    parser.add_argument("--genome", required=True, help="Genome FASTA path.")
    parser.add_argument("--output", required=True, help="Output CDS FASTA path.")
    parser.add_argument(
        "--blast-table",
        default=None,
        help="Optional output path for BLAST tabular results.",
    )
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold.")
    parser.add_argument("--threads", type=int, default=4, help="BLAST threads.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    proteins = Path(args.proteins)
    genome = Path(args.genome)
    output = Path(args.output)

    try:
        require_tool("makeblastdb")
        require_tool("tblastn")
    except RuntimeError as exc:
        logging.error(str(exc))
        return 1

    db_prefix = genome.with_suffix("")
    try:
        make_blast_db(genome, db_prefix)
    except subprocess.CalledProcessError as exc:
        logging.error("makeblastdb failed: %s", exc)
        return 1

    out_table = Path(args.blast_table) if args.blast_table else output.with_suffix(".tblastn.tsv")
    try:
        run_tblastn(proteins, db_prefix, out_table, args.evalue, args.threads)
    except subprocess.CalledProcessError as exc:
        logging.error("tblastn failed: %s", exc)
        return 1

    hits = parse_hits(out_table, args.evalue)
    if not hits:
        logging.error("No BLAST hits found with e-value <= %s", args.evalue)
        return 1

    genome_records = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    extracted = extract_cds(hits, genome_records)
    if not extracted:
        logging.error("No CDS sequences extracted.")
        return 1

    output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(extracted, output, "fasta")
    logging.info("Wrote %d CDS sequences to %s", len(extracted), output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
