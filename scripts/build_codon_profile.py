"""Build a FactorForge codon profile from CDS FASTA records."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_PATH = ROOT / "src" / "factorforge" / "data" / "nbenthamiana_codons.json"
STOP_CODONS = {"TAA", "TAG", "TGA"}
VALID_BASES = set("ATGC")
ORGANELLAR_MARKERS = ("chloroplast", "plastid", "mitochondria", "chrc", "chrm")
AA_NAMES = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic acid",
    "E": "Glutamic acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "*": "Stop",
}
STANDARD_GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str
    order: int


def _read_fasta(path: Path) -> list[FastaRecord]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(chunks), len(records)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append(FastaRecord(header, "".join(chunks), len(records)))
    return records


def _trim_terminal_stop(sequence: str) -> str:
    terminal = sequence[-3:] if len(sequence) >= 3 else ""
    if terminal in STOP_CODONS:
        return sequence[:-3]
    return sequence


def _has_internal_stop(sequence: str) -> bool:
    return any(sequence[index : index + 3] in STOP_CODONS for index in range(0, len(sequence), 3))


def _is_organellar(header: str) -> bool:
    lowered = header.lower()
    return any(marker in lowered for marker in ORGANELLAR_MARKERS)


def _extract_gene_id(header: str) -> str | None:
    for pattern in (r"(?:^|[\s;|])gene=([^\s;|]+)", r"(?:^|[\s;|])gene:([^\s;|]+)"):
        match = re.search(pattern, header, flags=re.IGNORECASE)
        if match:
            return match.group(1)

    first_token = header.split()[0] if header.split() else ""
    if first_token.lower().startswith("transcript:"):
        transcript = first_token.split(":", 1)[1]
        if not transcript:
            return None
        parts = transcript.rsplit(".", 1)
        return parts[0] if len(parts) == 2 and parts[1].isdigit() else transcript
    parts = first_token.rsplit(".", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return first_token or None


def _filter_records(records: list[FastaRecord]) -> tuple[list[FastaRecord], dict[str, int]]:
    stats = {
        "records_input": len(records),
        "exclusion_non_atgc": 0,
        "exclusion_length_mod3": 0,
        "exclusion_too_short": 0,
        "exclusion_internal_stop": 0,
        "exclusion_organellar": 0,
        "exclusion_duplicate": 0,
        "exclusion_isoform_shorter": 0,
    }
    filtered: list[FastaRecord] = []
    seen_hashes: set[str] = set()

    for record in records:
        sequence = record.sequence.upper()
        if set(sequence) - VALID_BASES:
            stats["exclusion_non_atgc"] += 1
            continue
        sequence = _trim_terminal_stop(sequence)
        if len(sequence) % 3 != 0:
            stats["exclusion_length_mod3"] += 1
            continue
        if len(sequence) < 60:
            stats["exclusion_too_short"] += 1
            continue
        if _has_internal_stop(sequence):
            stats["exclusion_internal_stop"] += 1
            continue
        if _is_organellar(record.header):
            stats["exclusion_organellar"] += 1
            continue
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if digest in seen_hashes:
            stats["exclusion_duplicate"] += 1
            continue
        seen_hashes.add(digest)
        filtered.append(FastaRecord(record.header, sequence, record.order))

    no_gene_id: list[FastaRecord] = []
    best_by_gene: dict[str, FastaRecord] = {}
    for record in filtered:
        gene_id = _extract_gene_id(record.header)
        if gene_id is None:
            no_gene_id.append(record)
            continue
        current = best_by_gene.get(gene_id)
        if current is None:
            best_by_gene[gene_id] = record
        elif len(record.sequence) > len(current.sequence):
            best_by_gene[gene_id] = record
            stats["exclusion_isoform_shorter"] += 1
        else:
            stats["exclusion_isoform_shorter"] += 1

    kept = sorted([*no_gene_id, *best_by_gene.values()], key=lambda record: record.order)
    return kept, stats


def _count_codons(records: list[FastaRecord]) -> tuple[dict[str, int], int, int]:
    counts = {codon: 0 for codon in STANDARD_GENETIC_CODE}
    gc_count = 0
    total_nt = 0
    total_sense = 0
    for record in records:
        gc_count += record.sequence.count("G") + record.sequence.count("C")
        total_nt += len(record.sequence)
        for index in range(0, len(record.sequence), 3):
            codon = record.sequence[index : index + 3]
            aa = STANDARD_GENETIC_CODE.get(codon)
            if aa is None:
                continue
            counts[codon] += 1
            if aa != "*":
                total_sense += 1
    return counts, total_sense, gc_count if total_nt else 0


def _group_codons_by_aa(codon_order: list[str]) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    for codon in codon_order:
        aa = STANDARD_GENETIC_CODE[codon]
        grouped.setdefault(aa, []).append(codon)
    return grouped


def _build_profile(
    source_profile_id: str,
    source_url: str,
    counts: dict[str, int],
    total_sense_codons: int,
    gc_count: int,
    total_nt: int,
) -> dict:
    template = json.loads(TEMPLATE_PATH.read_text(encoding="utf-8"))
    codon_order = list(template["codons"])
    codons_by_aa = _group_codons_by_aa(codon_order)
    codons = {}
    for codon in codon_order:
        aa = STANDARD_GENETIC_CODE[codon]
        if aa == "*":
            frequency = 0.0
            per_thousand = 0.0
        else:
            aa_codons = codons_by_aa[aa]
            denominator = sum(counts[item] for item in aa_codons) + len(aa_codons)
            frequency = (counts[codon] + 1) / denominator if denominator else 0.0
            per_thousand = counts[codon] / total_sense_codons * 1000 if total_sense_codons else 0.0
        codons[codon] = {
            "aa": aa,
            "frequency": round(frequency, 6),
            "per_thousand": round(per_thousand, 4),
        }

    amino_acids = {}
    for aa, metadata in template["amino_acids"].items():
        aa_codons = metadata["codons"]
        if aa == "*":
            preferred = metadata.get("preferred", "TAA")
        else:
            preferred = max(aa_codons, key=lambda codon: (codons[codon]["frequency"], counts[codon]))
        amino_acids[aa] = {
            "name": metadata.get("name", AA_NAMES[aa]),
            "codons": aa_codons,
            "preferred": preferred,
        }

    gc_fraction = gc_count / total_nt if total_nt else 0.0
    return {
        "organism": "Nicotiana benthamiana",
        "source": source_url,
        "description": (
            f"{source_profile_id} codon usage table built with strict_nuclear_cds_v1 "
            "filtering, Laplace pseudocount, build_codon_profile.py"
        ),
        "codons": codons,
        "amino_acids": amino_acids,
        "gc_content": {
            "overall": round(gc_fraction, 4),
            "description": "GC fraction across retained strict_nuclear_cds_v1 CDS records",
        },
        "notes": [
            (
                "Laplace pseudocount (count+1) applied. Filtering: strict_nuclear_cds_v1. "
                f"Build date: {date.today().isoformat()}."
            )
        ],
    }


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=4) + "\n", encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--source-profile-id", required=True)
    parser.add_argument("--source-url", required=True)
    parser.add_argument("--source-file-name", required=True)
    parser.add_argument("--download-date", required=True)
    parser.add_argument("--filtering-policy", required=True)
    parser.add_argument("--genetic-code", default="standard")
    parser.add_argument("--pseudocount-policy", required=True)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-manifest", required=True, type=Path)
    parser.add_argument("--out-filtered-stats", required=True, type=Path)
    args = parser.parse_args()
    if args.filtering_policy != "strict_nuclear_cds_v1":
        raise SystemExit("--filtering-policy must be strict_nuclear_cds_v1")
    if args.pseudocount_policy != "laplace":
        raise SystemExit("--pseudocount-policy must be laplace")
    if args.genetic_code != "standard":
        raise SystemExit("--genetic-code must be standard")
    if not args.input_fasta.exists():
        raise SystemExit(f"Input FASTA not found: {args.input_fasta}")
    return args


def main() -> None:
    args = _parse_args()
    records = _read_fasta(args.input_fasta)
    kept_records, filter_stats = _filter_records(records)
    counts, total_sense_codons, gc_count = _count_codons(kept_records)
    total_nt = sum(len(record.sequence) for record in kept_records)

    profile = _build_profile(
        args.source_profile_id,
        args.source_url,
        counts,
        total_sense_codons,
        gc_count,
        total_nt,
    )
    _write_json(args.out_json, profile)
    profile_sha256 = hashlib.sha256(args.out_json.read_bytes()).hexdigest()

    manifest = {
        "source_profile_id": args.source_profile_id,
        "source_url": args.source_url,
        "source_file_name": args.source_file_name,
        "download_date": args.download_date,
        "build_date": date.today().isoformat(),
        "filtering_policy": args.filtering_policy,
        "pseudocount_policy": args.pseudocount_policy,
        "genetic_code": args.genetic_code,
        "records_input": filter_stats["records_input"],
        "records_used": len(kept_records),
        "exclusion_non_atgc": filter_stats["exclusion_non_atgc"],
        "exclusion_length_mod3": filter_stats["exclusion_length_mod3"],
        "exclusion_too_short": filter_stats["exclusion_too_short"],
        "exclusion_internal_stop": filter_stats["exclusion_internal_stop"],
        "exclusion_organellar": filter_stats["exclusion_organellar"],
        "exclusion_duplicate": filter_stats["exclusion_duplicate"],
        "exclusion_isoform_shorter": filter_stats["exclusion_isoform_shorter"],
        "total_sense_codons_counted": total_sense_codons,
        "codon_profile_sha256": profile_sha256,
        "build_script": "scripts/build_codon_profile.py",
    }
    stats = {
        "source_profile_id": args.source_profile_id,
        "total_sense_codons_counted": total_sense_codons,
        "raw_codon_counts_before_pseudocount": counts,
    }
    _write_json(args.out_manifest, manifest)
    _write_json(args.out_filtered_stats, stats)

    print(
        f"{args.source_profile_id}: records_input={filter_stats['records_input']} "
        f"records_used={len(kept_records)} total_sense_codons={total_sense_codons} "
        f"sha256={profile_sha256}"
    )


if __name__ == "__main__":
    main()
