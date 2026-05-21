"""Compute codon identity between v3 decoded CDS and v2 gc_target candidates."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2  # noqa: E402


def _normalize_sequence(sequence: str) -> str:
    return "".join(sequence.upper().replace("U", "T").split())


def _codons(sequence: str) -> list[str]:
    seq = _normalize_sequence(sequence)
    end = len(seq) - len(seq) % 3
    return [seq[index : index + 3] for index in range(0, end, 3)]


def codon_identity_stats(v3_cds: str, v2_cds: str) -> tuple[float, int]:
    """Return percent identical codons and number of different codon positions."""
    v3_codons = _codons(v3_cds)
    v2_codons = _codons(v2_cds)
    if not v3_codons or len(v3_codons) != len(v2_codons):
        return 0.0, max(len(v3_codons), len(v2_codons))
    different_positions = sum(
        1 for v3_codon, v2_codon in zip(v3_codons, v2_codons) if v3_codon != v2_codon
    )
    identity_pct = ((len(v3_codons) - different_positions) / len(v3_codons)) * 100.0
    return identity_pct, different_positions


def _protein_sequence(row: dict[str, str]) -> str:
    sequence = row.get("amino_acid_sequence") or row.get("protein_sequence") or row.get("sequence")
    if not sequence:
        raise ValueError(f"Missing amino acid sequence for protein_id={row.get('protein_id')}")
    return "".join(sequence.upper().split()).rstrip("*")


def _v3_cds(row: dict[str, str]) -> str:
    sequence = row.get("v3_decoded_cds") or row.get("v3_cds") or row.get("v3_dna_sequence")
    if not sequence:
        raise ValueError(f"Missing v3 decoded CDS for protein_id={row.get('protein_id')}")
    return _normalize_sequence(sequence)


def read_v3_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_identity_rows(v3_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    identity_rows: list[dict[str, Any]] = []
    for row in v3_rows:
        protein_id = str(row["protein_id"])
        protein = _protein_sequence(row)
        v3_cds = _v3_cds(row)
        v2_result = optimize_with_v2(protein, options={"profile": "gc_target", "scan_mode": "fast"})
        v2_cds = _normalize_sequence(str(v2_result["dna_sequence"]))
        identity_pct, different_positions = codon_identity_stats(v3_cds, v2_cds)
        identity_rows.append(
            {
                "protein_id": protein_id,
                "v3_cds_length": len(v3_cds),
                "v2_gc_target_cds": v2_cds,
                "v3_cds": v3_cds,
                "codon_identity_pct": identity_pct,
                "different_positions": different_positions,
            }
        )
    return identity_rows


def summarize_identity(identity_rows: list[dict[str, Any]]) -> dict[str, float | int]:
    values = [float(row["codon_identity_pct"]) for row in identity_rows]
    if not values:
        raise ValueError("No identity rows to summarize")
    return {
        "count": len(values),
        "mean_codon_identity": statistics.fmean(values),
        "std": statistics.pstdev(values),
        "min": min(values),
        "max": max(values),
        "pct_above_90": sum(1 for value in values if value >= 90.0) / len(values) * 100.0,
        "pct_above_95": sum(1 for value in values if value >= 95.0) / len(values) * 100.0,
        "pct_above_99": sum(1 for value in values if value >= 99.0) / len(values) * 100.0,
    }


def write_identity_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "protein_id",
        "v3_cds_length",
        "v2_gc_target_cds",
        "v3_cds",
        "codon_identity_pct",
        "different_positions",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_identity(
    v3_decoded_csv: Path,
    output_csv: Path,
    summary_json: Path,
) -> dict[str, float | int]:
    rows = build_identity_rows(read_v3_rows(v3_decoded_csv))
    summary = summarize_identity(rows)
    write_identity_csv(output_csv, rows)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute v3 vs v2 gc_target codon identity.")
    parser.add_argument(
        "--v3-decoded",
        type=Path,
        default=ROOT / "experiments" / "results" / "alpha_run2" / "v3_decoded_eval.csv",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=ROOT
        / "experiments"
        / "results"
        / "alpha_run2"
        / "codon_identity_report.csv",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=ROOT
        / "experiments"
        / "results"
        / "alpha_run2"
        / "codon_identity_summary.json",
    )
    args = parser.parse_args()
    print(
        json.dumps(
            run_identity(args.v3_decoded, args.output_csv, args.summary_json),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
