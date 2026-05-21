#!/usr/bin/env python3
"""
Run Salmon quantification and merge quant.sf files.
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found on PATH: {tool}")


def run_command(cmd: List[str]) -> None:
    logging.debug("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def build_index(transcripts: Path, index_dir: Path) -> None:
    if index_dir.exists() and any(index_dir.iterdir()):
        logging.info("Salmon index exists: %s", index_dir)
        return
    index_dir.mkdir(parents=True, exist_ok=True)
    run_command(["salmon", "index", "-t", str(transcripts), "-i", str(index_dir)])


def discover_reads(reads_dir: Path) -> Dict[str, List[Path]]:
    reads: Dict[str, List[Path]] = {}
    for path in reads_dir.rglob("*"):
        is_fastq = path.is_file() and (
            path.suffix.lower() in {".fastq", ".fq"} or str(path).endswith((".fastq.gz", ".fq.gz"))
        )
        if is_fastq:
            name = path.name
            for suffix in FASTQ_SUFFIXES:
                if name.endswith(suffix):
                    sample = name[: -len(suffix)]
                    break
            else:
                sample = name
            sample = sample.replace("_1", "").replace("_2", "")
            reads.setdefault(sample, []).append(path)
    return reads


def find_pairs(paths: List[Path]) -> Tuple[List[Path], List[Path], List[Path]]:
    r1 = [p for p in paths if "_1" in p.name or ".1" in p.name]
    r2 = [p for p in paths if "_2" in p.name or ".2" in p.name]
    singles = [p for p in paths if p not in r1 and p not in r2]
    return r1, r2, singles


def run_quant(
    sample: str,
    reads: List[Path],
    index_dir: Path,
    output_dir: Path,
    libtype: str,
    threads: int,
) -> None:
    r1, r2, singles = find_pairs(reads)
    sample_out = output_dir / sample
    sample_out.mkdir(parents=True, exist_ok=True)

    if r1 and r2:
        run_command(
            [
                "salmon",
                "quant",
                "-i",
                str(index_dir),
                "-l",
                libtype,
                "-1",
                str(sorted(r1)[0]),
                "-2",
                str(sorted(r2)[0]),
                "-p",
                str(threads),
                "-o",
                str(sample_out),
            ]
        )
    elif singles:
        run_command(
            [
                "salmon",
                "quant",
                "-i",
                str(index_dir),
                "-l",
                libtype,
                "-r",
                str(sorted(singles)[0]),
                "-p",
                str(threads),
                "-o",
                str(sample_out),
            ]
        )
    else:
        raise RuntimeError(f"No readable FASTQ files for sample {sample}")


def merge_quants(output_dir: Path, merged_path: Path, metric: str) -> None:
    samples = []
    quant_files = []
    for sample_dir in sorted(output_dir.iterdir()):
        quant_file = sample_dir / "quant.sf"
        if sample_dir.is_dir() and quant_file.exists():
            samples.append(sample_dir.name)
            quant_files.append(quant_file)

    if not quant_files:
        raise RuntimeError("No quant.sf files found for merging.")

    merged: Dict[str, Dict[str, str]] = {}
    for sample, quant_file in zip(samples, quant_files):
        with quant_file.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if metric not in reader.fieldnames:
                raise RuntimeError(f"Metric {metric} not found in {quant_file}")
            for row in reader:
                name = row["Name"]
                merged.setdefault(name, {})[sample] = row[metric]

    merged_path.parent.mkdir(parents=True, exist_ok=True)
    with merged_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Name"] + samples)
        for name in sorted(merged.keys()):
            row = [name] + [merged[name].get(sample, "0") for sample in samples]
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Salmon quant and merge quant.sf files.")
    parser.add_argument("--transcripts", required=True, help="Transcriptome FASTA path.")
    parser.add_argument("--reads-dir", required=True, help="Directory with FASTQ files.")
    parser.add_argument("--index", required=True, help="Salmon index directory.")
    parser.add_argument(
        "--output",
        default="data/rnaseq/salmon_quant",
        help="Output directory for Salmon quant results.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Salmon threads.")
    parser.add_argument("--libtype", default="A", help="Salmon library type.")
    parser.add_argument("--metric", default="TPM", help="Metric to merge from quant.sf.")
    parser.add_argument(
        "--merged-output",
        default="data/rnaseq/salmon_quant/merged_quant.tsv",
        help="Merged quant output path.",
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    transcripts = Path(args.transcripts)
    reads_dir = Path(args.reads_dir)
    index_dir = Path(args.index)
    output_dir = Path(args.output)
    merged_output = Path(args.merged_output)

    try:
        require_tool("salmon")
    except RuntimeError as exc:
        logging.error(str(exc))
        return 1

    try:
        build_index(transcripts, index_dir)
    except subprocess.CalledProcessError as exc:
        logging.error("Salmon index failed: %s", exc)
        return 1

    reads = discover_reads(reads_dir)
    if not reads:
        logging.error("No FASTQ files found under %s", reads_dir)
        return 1

    for sample, paths in sorted(reads.items()):
        try:
            logging.info("Quantifying sample %s", sample)
            run_quant(sample, paths, index_dir, output_dir, args.libtype, args.threads)
        except Exception as exc:
            logging.error("Quant failed for %s: %s", sample, exc)
            return 1

    try:
        merge_quants(output_dir, merged_output, args.metric)
    except Exception as exc:
        logging.error("Merge failed: %s", exc)
        return 1

    logging.info("Merged quant written to %s", merged_output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
