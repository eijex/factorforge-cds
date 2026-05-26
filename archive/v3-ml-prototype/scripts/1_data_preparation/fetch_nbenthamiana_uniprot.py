"""Fetch N. benthamiana proteins from the UniProt REST API."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"


def _next_link(link_header: str) -> str | None:
    for part in link_header.split(","):
        if 'rel="next"' in part and "<" in part and ">" in part:
            return part.split("<", 1)[1].split(">", 1)[0]
    return None


def fetch_sequences(
    max_sequences: int = 15000,
    output_path: str = "data/raw/uniprot_nbenthamiana_extended.fasta",
) -> None:
    """Download up to ``max_sequences`` unreviewed N. benthamiana proteins."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    params = {
        "query": "organism_id:4100 AND reviewed:false",
        "format": "fasta",
        "size": 500,
    }

    total = 0
    next_url: str | None = UNIPROT_API
    with out.open("w", encoding="utf-8") as handle:
        while next_url and total < max_sequences:
            if next_url == UNIPROT_API:
                response = requests.get(next_url, params=params, timeout=30)
            else:
                response = requests.get(next_url, timeout=30)
            response.raise_for_status()

            entries = [entry for entry in response.text.split(">") if entry.strip()]
            for entry in entries:
                if total >= max_sequences:
                    break
                handle.write(">" + entry)
                total += 1

            next_url = _next_link(response.headers.get("Link", ""))
            print(f"Downloaded: {total}", end="\r")
            time.sleep(0.5)

    print(f"\nSaved {total} sequences to {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max", type=int, default=15000)
    parser.add_argument("--output", default="data/raw/uniprot_nbenthamiana_extended.fasta")
    args = parser.parse_args()
    fetch_sequences(args.max, args.output)


if __name__ == "__main__":
    main()
