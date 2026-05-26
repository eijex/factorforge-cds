"""Fetch and filter N. benthamiana CDS FASTA files from the SGN FTP index."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BASE_URL = "https://solgenomics.net/ftp/genomes/Nicotiana_benthamiana/"
DEFAULT_OUTPUT = Path("data/raw/sgn_nbenthamiana_cds.fasta")
STOP_CODONS = {"TAA", "TAG", "TGA"}


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.links.append(value)


def parse_links(html: str) -> list[str]:
    parser = LinkParser()
    parser.feed(html)
    return parser.links


def fetch_index(url: str) -> list[str]:
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    return parse_links(response.text)


def find_cds_fasta_url(base_url: str = BASE_URL) -> str:
    """Find the best SGN CDS FASTA URL under the N. benthamiana FTP index."""
    candidate_urls: list[str] = []
    queue = [base_url]
    seen: set[str] = set()

    while queue:
        url = queue.pop(0)
        if url in seen:
            continue
        seen.add(url)
        for href in fetch_index(url):
            if href in {"../", "/"}:
                continue
            child_url = urljoin(url, href)
            name = href.lower()
            if href.endswith("/") and any(token in child_url for token in ("LAB360", "QLD183", "annotation")):
                queue.append(child_url)
                continue
            if re.search(r"(\.cds\.fa|\.cds\.fasta|cds\.fasta)(\.gz)?$", name) and ".aa." not in name:
                candidate_urls.append(child_url)

    if not candidate_urls:
        raise RuntimeError(f"No CDS FASTA file found under {base_url}")

    # Prefer LAB360 when present; it is the current SGN N. benthamiana line used by this job.
    return sorted(candidate_urls, key=lambda item: ("LAB360" not in item, item))[0]


def iter_fasta_text(url: str):
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()
    response.raw.decode_content = True
    if url.endswith(".gz"):
        with gzip.GzipFile(fileobj=response.raw) as handle:
            for raw_line in handle:
                yield raw_line.decode("utf-8")
    else:
        for raw_line in response.iter_lines(decode_unicode=True):
            yield raw_line + "\n"


def parse_fasta_lines(lines) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    current_id: str | None = None
    current_description = ""
    current_lines: list[str] = []

    def flush() -> None:
        if current_id is None:
            return
        records.append(
            SeqRecord(
                Seq("".join(current_lines)),
                id=current_id,
                description=current_description,
            )
        )

    for line in lines:
        text = line.strip()
        if not text:
            continue
        if text.startswith(">"):
            flush()
            current_description = text[1:]
            current_id = current_description.split()[0]
            current_lines = []
        else:
            current_lines.append(text)
    flush()
    return records


def is_valid_orf(seq: str, min_len: int = 150, max_len: int = 5000) -> bool:
    seq = seq.upper().replace("U", "T")
    if len(seq) < min_len or len(seq) > max_len or len(seq) % 3 != 0:
        return False
    if set(seq) - set("ATGC"):
        return False
    return seq.startswith("ATG") and seq[-3:] in STOP_CODONS


def fetch_records(url: str | None = None) -> list[SeqRecord]:
    cds_url = url or find_cds_fasta_url()
    records: list[SeqRecord] = []
    for record in parse_fasta_lines(iter_fasta_text(cds_url)):
        seq = str(record.seq).upper().replace("U", "T")
        if is_valid_orf(seq):
            record.seq = Seq(seq)
            records.append(record)
    return records


def write_records(records: list[SeqRecord], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output, "fasta")


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch SGN N. benthamiana CDS FASTA.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--url", default=None, help="Override the auto-discovered CDS FASTA URL.")
    parser.add_argument("--count", action="store_true", help="Print filtered CDS count only.")
    args = parser.parse_args()

    try:
        cds_url = args.url or find_cds_fasta_url()
        records = fetch_records(cds_url)
    except requests.RequestException as exc:
        print(f"FTP connection required or unavailable: {exc}")
        return 0 if args.count else 1
    except RuntimeError as exc:
        print(f"FTP connection required or unavailable: {exc}")
        return 0 if args.count else 1

    if args.count:
        print(len(records))
        return 0

    write_records(records, args.output)
    print(f"Source: {cds_url}")
    print(f"Saved {len(records)} filtered CDS records to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
