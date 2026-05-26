#!/usr/bin/env python3
"""
Download RNA-seq data from NCBI SRA using SRA Toolkit.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List

DEFAULT_ACCESSIONS = ["SRR5983857", "SRR5983858", "SRR5983859", "SRR5983860"]


def setup_logging(verbose: bool, log_file: Path) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(levelname)s: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode="a", encoding="utf-8"),
        ],
    )


def require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found on PATH: {tool}")


def run_command(cmd: List[str]) -> None:
    logging.debug("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def find_sra_file(sra_dir: Path, accession: str) -> Path:
    candidate_dir = sra_dir / accession
    for path in candidate_dir.rglob("*.sra"):
        return path
    raise FileNotFoundError(f"SRA file not found for {accession} in {candidate_dir}")


def download_accession(accession: str, sra_dir: Path, fastq_dir: Path, threads: int) -> None:
    run_command(["prefetch", accession, "--output-directory", str(sra_dir)])
    sra_file = find_sra_file(sra_dir, accession)
    run_command(
        [
            "fasterq-dump",
            str(sra_file),
            "--split-files",
            "--threads",
            str(threads),
            "--outdir",
            str(fastq_dir),
        ]
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download SRP109347 RNA-seq data.")
    parser.add_argument("--project", default="SRP109347", help="SRA project ID.")
    parser.add_argument(
        "--accessions",
        nargs="*",
        default=DEFAULT_ACCESSIONS,
        help="SRA run accessions.",
    )
    parser.add_argument(
        "--output",
        default="data/rnaseq/SRP109347",
        help="Output directory for SRA and FASTQ files.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Threads for fasterq-dump.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    base_dir = Path(args.output)
    sra_dir = base_dir / "sra"
    fastq_dir = base_dir / "fastq"
    base_dir.mkdir(parents=True, exist_ok=True)

    log_file = base_dir / "download_log.txt"
    setup_logging(args.verbose, log_file)

    try:
        require_tool("prefetch")
        require_tool("fasterq-dump")
    except RuntimeError as exc:
        logging.error(str(exc))
        return 1

    logging.info("Project: %s", args.project)
    logging.info("Accessions: %s", ", ".join(args.accessions))

    for accession in args.accessions:
        try:
            logging.info("Downloading %s", accession)
            download_accession(accession, sra_dir, fastq_dir, args.threads)
        except Exception as exc:
            logging.error("Failed %s: %s", accession, exc)
            return 1

    logging.info("Download complete: %s", base_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
