#!/usr/bin/env python3
"""
Run CD-HIT-EST to reduce redundancy and summarize cluster statistics.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import mean
from typing import List


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found on PATH: {tool}")


def run_command(cmd: List[str]) -> None:
    logging.debug("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def choose_word_size(identity: float) -> int:
    if identity >= 0.95:
        return 10
    if identity >= 0.9:
        return 8
    if identity >= 0.88:
        return 7
    return 6


def parse_clusters(cluster_path: Path) -> List[int]:
    clusters: List[int] = []
    current = 0
    with cluster_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = 0
            else:
                current += 1
    if current:
        clusters.append(current)
    return clusters


def write_stats(stats_path: Path, clusters: List[int]) -> None:
    total_clusters = len(clusters)
    total_sequences = sum(clusters)
    max_size = max(clusters) if clusters else 0
    min_size = min(clusters) if clusters else 0
    avg_size = mean(clusters) if clusters else 0.0
    singletons = sum(1 for c in clusters if c == 1)

    stats_path.parent.mkdir(parents=True, exist_ok=True)
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write("CD-HIT-EST CLUSTER SUMMARY\n")
        handle.write(f"Total sequences: {total_sequences}\n")
        handle.write(f"Total clusters: {total_clusters}\n")
        handle.write(f"Singleton clusters: {singletons}\n")
        handle.write(f"Min cluster size: {min_size}\n")
        handle.write(f"Max cluster size: {max_size}\n")
        handle.write(f"Average cluster size: {avg_size:.2f}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CD-HIT-EST and summarize clusters.")
    parser.add_argument("--input", required=True, help="Input FASTA path.")
    parser.add_argument("--output", required=True, help="Output FASTA path.")
    parser.add_argument("--identity", type=float, default=0.95, help="Identity threshold.")
    parser.add_argument("--threads", type=int, default=4, help="Threads for CD-HIT-EST.")
    parser.add_argument("--memory", type=int, default=0, help="Memory limit (MB). 0=unlimited.")
    parser.add_argument("--word-size", type=int, default=None, help="Word size for CD-HIT-EST.")
    parser.add_argument(
        "--stats",
        default=None,
        help="Output path for cluster statistics.",
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    input_path = Path(args.input)
    output_path = Path(args.output)
    stats_path = Path(args.stats) if args.stats else output_path.with_suffix(".stats.txt")
    word_size = args.word_size or choose_word_size(args.identity)

    try:
        require_tool("cd-hit-est")
    except RuntimeError as exc:
        logging.error(str(exc))
        return 1

    cmd = [
        "cd-hit-est",
        "-i",
        str(input_path),
        "-o",
        str(output_path),
        "-c",
        str(args.identity),
        "-n",
        str(word_size),
        "-T",
        str(args.threads),
        "-M",
        str(args.memory),
    ]

    try:
        run_command(cmd)
    except subprocess.CalledProcessError as exc:
        logging.error("CD-HIT-EST failed: %s", exc)
        return 1

    cluster_path = output_path.with_suffix(".clstr")
    if not cluster_path.exists():
        logging.error("Cluster file not found: %s", cluster_path)
        return 1

    clusters = parse_clusters(cluster_path)
    write_stats(stats_path, clusters)
    logging.info("Cluster stats written to %s", stats_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
