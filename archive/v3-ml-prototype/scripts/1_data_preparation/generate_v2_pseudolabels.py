"""Generate v2 pseudo-labels for v3-alpha Run 4.

This script is intentionally data-preparation only. It does not train v3.
Each output row carries the protein sequence, v2 pseudo-label CDS, spaced codons,
metrics, and structured validator output so unsafe rows can be filtered before
training.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path
from typing import Any, Iterable

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2


def _spaced_codons(dna_sequence: str) -> str:
    seq = "".join(dna_sequence.upper().split())
    return " ".join(seq[index : index + 3] for index in range(0, len(seq), 3))


def _protein_records_from_fasta(fasta_path: Path) -> Iterable[tuple[str, str]]:
    for record in SeqIO.parse(fasta_path, "fasta"):
        yield record.id, str(record.seq).upper().replace("*", "")


def _protein_records_from_cds_fasta(fasta_path: Path) -> Iterable[tuple[str, str]]:
    for record in SeqIO.parse(fasta_path, "fasta"):
        dna = str(record.seq).upper().replace("U", "T").replace(" ", "")
        if len(dna) < 3 or len(dna) % 3 != 0 or set(dna) - set("ATGC"):
            continue
        protein = str(Seq(dna).translate(to_stop=False)).rstrip("*")
        if "*" in protein:
            continue
        if protein:
            yield record.id, protein


def _write_row(handle, row: dict[str, Any]) -> None:
    handle.write(json.dumps(row, sort_keys=True) + "\n")


def _build_row(
    protein_id: str,
    protein_sequence: str,
    profile: str,
    require_validator_pass: bool,
) -> dict[str, Any] | None:
    result = optimize_with_v2(
        protein_sequence,
        options={"profile": profile, "scan_mode": "fast"},
    )
    validator = result["validator"]
    if require_validator_pass and not validator["passed"]:
        return None

    return {
        "protein_id": protein_id,
        "sequence": result["protein_sequence"],
        "protein_sequence": result["protein_sequence"],
        "amino_acid_sequence": result["protein_sequence"],
        "codon_sequence": _spaced_codons(result["dna_sequence"]),
        "dna_seq": result["dna_sequence"],
        "dna_sequence": result["dna_sequence"],
        "pseudo_label_engine": "v2",
        "pseudo_label_profile": profile,
        "metrics": result["metrics"],
        "validator": validator,
    }


def generate_v2_pseudolabels(
    fasta_path: str | Path,
    output_path: str | Path,
    profile: str = "high_cai",
    max_length: int = 512,
    require_validator_pass: bool = True,
    skip_on_error: bool = True,
    input_kind: str = "protein",
    limit: int | None = None,
) -> dict[str, int]:
    """Generate v2 pseudo-label JSONL rows from a protein FASTA file."""
    fasta = Path(fasta_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    counts = {
        "read": 0,
        "written": 0,
        "skipped_long": 0,
        "skipped_invalid": 0,
        "skipped_error": 0,
    }

    if input_kind not in {"protein", "cds"}:
        raise ValueError("input_kind must be 'protein' or 'cds'")
    records = (
        _protein_records_from_cds_fasta(fasta)
        if input_kind == "cds"
        else _protein_records_from_fasta(fasta)
    )

    with output.open("w", encoding="utf-8") as handle:
        for protein_id, protein_sequence in records:
            counts["read"] += 1
            if limit is not None and counts["written"] >= limit:
                break
            protein = "".join(protein_sequence.upper().split()).rstrip("*")
            if len(protein) > max_length:
                protein = protein[:max_length]
                counts["skipped_long"] += 1
            try:
                row = _build_row(
                    protein_id=protein_id,
                    protein_sequence=protein,
                    profile=profile,
                    require_validator_pass=require_validator_pass,
                )
                if row is None:
                    counts["skipped_invalid"] += 1
                    continue
                _write_row(handle, row)
                counts["written"] += 1
            except Exception:
                counts["skipped_error"] += 1
                if not skip_on_error:
                    raise

    return counts


def split_jsonl(
    input_path: str | Path,
    train_output: str | Path,
    eval_output: str | Path,
    train_split: float = 0.9,
    seed: int = 42,
) -> dict[str, int]:
    """Split a pseudo-label JSONL file into deterministic train/eval files."""
    source = Path(input_path)
    rows = [line for line in source.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not 0.0 < train_split < 1.0:
        raise ValueError("train_split must be between 0 and 1")
    rng = random.Random(seed)
    rng.shuffle(rows)
    train_count = int(len(rows) * train_split)
    train_rows = rows[:train_count]
    eval_rows = rows[train_count:]

    train_path = Path(train_output)
    eval_path = Path(eval_output)
    train_path.parent.mkdir(parents=True, exist_ok=True)
    eval_path.parent.mkdir(parents=True, exist_ok=True)
    train_path.write_text("\n".join(train_rows) + ("\n" if train_rows else ""), encoding="utf-8")
    eval_path.write_text("\n".join(eval_rows) + ("\n" if eval_rows else ""), encoding="utf-8")
    return {"total": len(rows), "train": len(train_rows), "eval": len(eval_rows)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Run 4 v2 pseudo-label JSONL.")
    parser.add_argument("--fasta", required=True, help="Input protein FASTA.")
    parser.add_argument("--output", default="data/training/run4_v2_pseudolabels.jsonl")
    parser.add_argument("--train-output", default=None, help="Optional train split JSONL output.")
    parser.add_argument("--eval-output", default=None, help="Optional eval split JSONL output.")
    parser.add_argument("--train-split", type=float, default=0.9)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--limit", type=int, default=None, help="Maximum number of written rows.")
    parser.add_argument(
        "--input-kind",
        choices=["protein", "cds"],
        default="protein",
        help="Whether --fasta contains protein or CDS records.",
    )
    parser.add_argument(
        "--profile",
        default="high_cai",
        choices=["balanced", "high_cai", "gc_target", "assembly_friendly", "ramp", "viral_delivery"],
    )
    parser.add_argument("--max-length", type=int, default=512)
    parser.add_argument("--allow-validator-fail", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    args = parser.parse_args()

    counts = generate_v2_pseudolabels(
        fasta_path=args.fasta,
        output_path=args.output,
        profile=args.profile,
        max_length=args.max_length,
        require_validator_pass=not args.allow_validator_fail,
        skip_on_error=not args.fail_fast,
        input_kind=args.input_kind,
        limit=args.limit,
    )
    payload: dict[str, Any] = {"generation": counts}
    if args.train_output or args.eval_output:
        if not args.train_output or not args.eval_output:
            parser.error("--train-output and --eval-output must be provided together")
        payload["split"] = split_jsonl(
            input_path=args.output,
            train_output=args.train_output,
            eval_output=args.eval_output,
            train_split=args.train_split,
            seed=args.seed,
        )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
