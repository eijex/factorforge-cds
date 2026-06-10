"""Fetch + filter the N. benthamiana reference CDS dataset.
Produces nbenthamiana_reference_cds.fasta + nbenthamiana_reference_proteins.fasta.
Run manually (network required); raw data is NOT committed (.gitignore).

Usage:
    python fetch_dataset.py --url <verified_url>   # download from URL
    python fetch_dataset.py --file <local_path>    # use already-downloaded file
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path

from factorforge.analysis.metrics import translate_dna

DATASETS_DIR = Path(__file__).resolve().parent
CDS_OUT = DATASETS_DIR / "nbenthamiana_reference_cds.fasta"
PROT_OUT = DATASETS_DIR / "nbenthamiana_reference_proteins.fasta"
# Real source URL + license MUST be verified before use; pass via --url.
SOURCE_URL = ""


def _read_fasta(text: str) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf).upper()
            name, buf = line[1:].strip().split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf).upper()
    return seqs


def _filter(raw: dict[str, str]) -> dict[str, str]:
    """Apply the filters documented in dataset_provenance.md."""
    kept, seen = {}, set()
    for name, seq in raw.items():
        if set(seq) - set("ACGT"):          # 1. ACGT only
            continue
        if len(seq) % 3 != 0:                # 2. length multiple of 3
            continue
        prot = translate_dna(seq)
        if "*" in prot[:-1]:                 # 3. no internal stop
            continue
        if not prot.rstrip("*"):             # 4. non-empty protein after stop removal
            continue
        h = hashlib.sha256(seq.encode()).hexdigest()
        if h in seen:                        # 4. dedupe by sequence hash
            continue
        seen.add(h)
        kept[name] = seq
    return kept


def _write_fasta(path: Path, seqs: dict[str, str]) -> str:
    import tempfile
    content = "".join(f">{n}\n{s}\n" for n, s in seqs.items()).encode("utf-8")
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False, suffix=".tmp") as tmp:
        tmp.write(content)
        tmp_path = Path(tmp.name)
    tmp_path.replace(path)  # atomic replace — prevents partial-file reads on failure
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=SOURCE_URL, required=False,
                    help="URL to a verified, redistributable N. benthamiana CDS FASTA.")
    ap.add_argument("--file", default=None, required=False,
                    help="Path to a locally downloaded CDS FASTA (plain text or .gz).")
    args = ap.parse_args()

    if args.file:
        p = Path(args.file)
        if not p.exists():
            raise SystemExit(f"File not found: {p}")
        if p.suffix == ".gz":
            import gzip
            raw_text = gzip.open(p, "rt", encoding="utf-8").read()
        else:
            raw_text = p.read_text(encoding="utf-8")
    elif args.url:
        import urllib.request
        resp = urllib.request.urlopen(args.url)  # noqa: S310
        if resp.status != 200:
            raise SystemExit(f"HTTP {resp.status} fetching dataset URL. Aborting.")
        raw_text = resp.read().decode("utf-8")
    else:
        raise SystemExit("Provide --url or --file. See module docstring for usage.")

    cds = _filter(_read_fasta(raw_text))                 # 5. CDS set
    proteins = {n: translate_dna(s).rstrip("*") for n, s in cds.items()}  # 6. protein set
    cds_hash = _write_fasta(CDS_OUT, cds)
    prot_hash = _write_fasta(PROT_OUT, proteins)
    (DATASETS_DIR / "checksums.txt").write_text(
        f"{CDS_OUT.name}  {cds_hash}\n{PROT_OUT.name}  {prot_hash}\n", encoding="utf-8")
    print(f"Kept {len(cds)} CDS. cds_sha256={cds_hash} prot_sha256={prot_hash}")


if __name__ == "__main__":
    main()
