#!/usr/bin/env python3
"""
Download N. benthamiana genome and annotation from Sol Genomics Network.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterable, Tuple

import requests
from tqdm import tqdm


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def normalize_organism(organism: str) -> str:
    org = organism.strip().lower()
    if org in {"nbenthamiana", "n.benthamiana", "nicotiana_benthamiana"}:
        return "Nicotiana_benthamiana"
    return organism


def candidate_urls(organism: str, version: str) -> Iterable[Tuple[str, str]]:
    org_folder = normalize_organism(organism)
    version_tag = version.lstrip("vV")
    bases = [
        f"https://ftp.solgenomics.net/genomes/{org_folder}/{version_tag}/",
        f"https://ftp.solgenomics.net/genomes/{org_folder}/",
    ]
    fasta_names = [
        f"{org_folder}_v{version_tag}_genome.fasta.gz",
        f"{org_folder}_v{version_tag}.fasta.gz",
        f"{org_folder}_{version_tag}.fasta.gz",
        f"{org_folder}.fasta.gz",
    ]
    gff_names = [
        f"{org_folder}_v{version_tag}_genome.gff3.gz",
        f"{org_folder}_v{version_tag}.gff3.gz",
        f"{org_folder}_{version_tag}.gff3.gz",
        f"{org_folder}.gff3.gz",
    ]

    for base in bases:
        for fasta_name in fasta_names:
            for gff_name in gff_names:
                yield base + fasta_name, base + gff_name


def download_with_progress(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Downloading %s", url)
    with requests.get(url, stream=True, timeout=60) as response:
        if response.status_code != 200:
            raise RuntimeError(f"Download failed with HTTP {response.status_code}: {url}")
        total = int(response.headers.get("Content-Length", 0))
        with (
            dest.open("wb") as handle,
            tqdm(total=total, unit="B", unit_scale=True, desc=dest.name) as bar,
        ):
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                handle.write(chunk)
                bar.update(len(chunk))


def resolve_urls(
    organism: str,
    version: str,
    fasta_url: str | None,
    gff_url: str | None,
) -> Tuple[str, str]:
    if fasta_url and gff_url:
        return fasta_url, gff_url
    if fasta_url or gff_url:
        raise ValueError("Provide both --fasta-url and --gff-url, or neither.")

    for fasta_candidate, gff_candidate in candidate_urls(organism, version):
        try:
            with requests.get(fasta_candidate, stream=True, timeout=15) as resp:
                if resp.status_code != 200:
                    continue
            with requests.get(gff_candidate, stream=True, timeout=15) as resp:
                if resp.status_code != 200:
                    continue
            return fasta_candidate, gff_candidate
        except requests.RequestException:
            continue

    raise RuntimeError(
        "Automatic URL resolution failed. Provide --fasta-url and --gff-url explicitly."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download N. benthamiana genome and GFF.")
    parser.add_argument("--organism", default="nbenthamiana", help="Organism name.")
    parser.add_argument("--version", default="v1.0.1", help="Genome version tag.")
    parser.add_argument(
        "--output",
        default="data/genomes/N_benthamiana_v1.0.1",
        help="Output directory.",
    )
    parser.add_argument("--fasta-url", default=None, help="Override FASTA URL.")
    parser.add_argument("--gff-url", default=None, help="Override GFF URL.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    output_dir = Path(args.output)
    fasta_url, gff_url = resolve_urls(args.organism, args.version, args.fasta_url, args.gff_url)

    fasta_name = Path(fasta_url).name
    gff_name = Path(gff_url).name

    try:
        download_with_progress(fasta_url, output_dir / fasta_name)
        download_with_progress(gff_url, output_dir / gff_name)
    except Exception as exc:
        logging.error(str(exc))
        return 1

    logging.info("Download complete: %s", output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
