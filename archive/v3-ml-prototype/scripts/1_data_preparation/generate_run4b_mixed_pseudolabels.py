"""Generate mixed Run 4B pseudo-labels for v3-alpha recovery.

Run 4 learned the low-GC high_cai teacher too closely. This generator builds
candidate teachers per protein and selects a validated multi-objective target
that can favor GC-feasible sequences without changing v2 production behavior.
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

from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile
from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2
from factorforge.engines.v3.metrics import load_codon_usage_table
from factorforge.ml.feasibility import analyze_feasibility
from factorforge.ml.metrics import (
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    detect_forbidden_motifs,
    detect_homopolymers,
    detect_repeats,
)
from factorforge.utils.validation import validate_candidate_sequence


GC_LOW = 40.0
GC_HIGH = 55.0
AVAILABLE_V2_PROFILES = {profile.value for profile in OptimizationProfile}


def _spaced_codons(dna_sequence: str) -> str:
    return " ".join(dna_sequence[index : index + 3] for index in range(0, len(dna_sequence), 3))


def _normalize_protein(sequence: str) -> str:
    return "".join(sequence.upper().split()).rstrip("*")


def _protein_records_from_fasta(fasta_path: Path) -> Iterable[tuple[str, str]]:
    for record in SeqIO.parse(fasta_path, "fasta"):
        protein = _normalize_protein(str(record.seq))
        if protein:
            yield record.id, protein


def _protein_records_from_cds_fasta(fasta_path: Path) -> Iterable[tuple[str, str]]:
    for record in SeqIO.parse(fasta_path, "fasta"):
        dna = str(record.seq).upper().replace("U", "T").replace(" ", "")
        if len(dna) < 3 or len(dna) % 3 != 0 or set(dna) - set("ATGC"):
            continue
        protein = str(Seq(dna).translate(to_stop=False)).rstrip("*")
        if protein and "*" not in protein:
            yield record.id, protein


def _protein_records_from_jsonl(jsonl_path: Path) -> Iterable[tuple[str, str]]:
    with jsonl_path.open(encoding="utf-8") as handle:
        for index, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            protein_id = str(row.get("protein_id") or f"row_{index:06d}")
            sequence = (
                row.get("amino_acid_sequence")
                or row.get("protein_sequence")
                or row.get("sequence")
            )
            if isinstance(sequence, str):
                protein = _normalize_protein(sequence)
                if protein:
                    yield protein_id, protein


def iter_protein_records(input_path: Path, input_kind: str) -> Iterable[tuple[str, str]]:
    if input_kind == "protein":
        return _protein_records_from_fasta(input_path)
    if input_kind == "cds":
        return _protein_records_from_cds_fasta(input_path)
    if input_kind == "jsonl":
        return _protein_records_from_jsonl(input_path)
    raise ValueError("input_kind must be 'protein', 'cds', or 'jsonl'")


def _metric_summary(protein: str, dna_sequence: str, codon_weights: dict[str, float]) -> dict[str, Any]:
    validator = validate_candidate_sequence(protein, dna_sequence)
    windows = calculate_gc_windows(dna_sequence)
    window_values = [float(window["gc"]) for window in windows]
    forbidden = detect_forbidden_motifs(dna_sequence, [])
    homopolymers = detect_homopolymers(dna_sequence)
    repeats = detect_repeats(dna_sequence)
    gc = calculate_gc(dna_sequence)
    cai = calculate_cai(dna_sequence, codon_weights)
    in_gc_range = GC_LOW <= gc <= GC_HIGH
    forbidden_count = len(forbidden)
    repeat_penalty_count = len(repeats) + len(homopolymers)
    selection_score = (
        cai
        + (0.20 if in_gc_range else 0.0)
        - (0.05 * forbidden_count)
        - (0.01 * repeat_penalty_count)
    )
    return {
        "dna_sequence": dna_sequence,
        "cai": cai,
        "gc": gc,
        "gc_in_40_55": in_gc_range,
        "gc_window_min": min(window_values) if window_values else 0.0,
        "gc_window_max": max(window_values) if window_values else 0.0,
        "forbidden_motif_count": forbidden_count,
        "homopolymer_count": len(homopolymers),
        "repeat_count": len(repeats),
        "selection_score": selection_score,
        "validator": validator,
    }


def _v2_candidate(
    protein: str,
    profile: str,
    codon_weights: dict[str, float],
    *,
    target_gc: float | None = None,
) -> dict[str, Any]:
    options: dict[str, Any] = {"profile": profile, "scan_mode": "fast"}
    if target_gc is not None:
        options["target_gc"] = target_gc
    result = optimize_with_v2(protein, options=options)
    metrics = _metric_summary(protein, result["dna_sequence"], codon_weights)
    return {"teacher_type": profile, **metrics}


def _feasibility_candidate(protein: str, codon_weights: dict[str, float]) -> dict[str, Any] | None:
    feasibility = analyze_feasibility(
        protein,
        codon_weights,
        target_gc_low=GC_LOW,
        target_gc_high=GC_HIGH,
        gc_ranges=[(GC_LOW, GC_HIGH)],
    )
    range_result = feasibility["ranges"].get(f"{GC_LOW:g}-{GC_HIGH:g}")
    if not range_result or not range_result.get("feasible"):
        return None
    best = range_result.get("best_candidate") or {}
    dna_sequence = best.get("dna_sequence")
    if not isinstance(dna_sequence, str) or not dna_sequence:
        return None
    metrics = _metric_summary(protein, dna_sequence, codon_weights)
    return {"teacher_type": "feasibility_40_55", **metrics}


def build_mixed_teacher_row(
    protein_id: str,
    protein_sequence: str,
    *,
    codon_weights: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Build one validated Run 4B mixed teacher row."""
    protein = _normalize_protein(protein_sequence)
    if not protein:
        raise ValueError("protein_sequence must not be empty")

    weights = codon_weights or load_codon_usage_table().codon_weights
    candidates: dict[str, dict[str, Any]] = {}

    high_cai = _v2_candidate(protein, "high_cai", weights)
    candidates["high_cai"] = high_cai

    if "gc_target" in AVAILABLE_V2_PROFILES:
        gc_target = _v2_candidate(protein, "gc_target", weights, target_gc=42.5)
        candidates["gc_target"] = gc_target
        if not gc_target["validator"]["passed"] or not gc_target["gc_in_40_55"]:
            feasibility = _feasibility_candidate(protein, weights)
            if feasibility is not None:
                candidates["feasibility_40_55"] = feasibility
    else:
        feasibility = _feasibility_candidate(protein, weights)
        if feasibility is not None:
            candidates["feasibility_40_55"] = feasibility

    valid_candidates = {
        name: candidate
        for name, candidate in candidates.items()
        if bool(candidate["validator"]["passed"])
    }
    if not valid_candidates:
        raise ValueError("No valid teacher candidates")

    selected_name, selected = max(
        valid_candidates.items(),
        key=lambda item: float(item[1]["selection_score"]),
    )
    dna_sequence = str(selected["dna_sequence"])
    validator = validate_candidate_sequence(protein, dna_sequence)
    candidate_metrics = {
        name: {
            key: value
            for key, value in candidate.items()
            if key not in {"dna_sequence"}
        }
        for name, candidate in candidates.items()
    }

    return {
        "protein_id": protein_id,
        "sequence": protein,
        "protein_sequence": protein,
        "amino_acid_sequence": protein,
        "codon_sequence": _spaced_codons(dna_sequence),
        "dna_seq": dna_sequence,
        "dna_sequence": dna_sequence,
        "selected_teacher_type": selected_name,
        "candidate_metrics": candidate_metrics,
        "selection_score": float(selected["selection_score"]),
        "validator_result": validator,
        "validator": validator,
        "pseudo_label_engine": "run4b_mixed",
        "pseudo_label_profile": selected_name,
    }


def generate_run4b_mixed_pseudolabels(
    input_path: str | Path,
    output_path: str | Path,
    *,
    input_kind: str = "protein",
    max_length: int = 512,
    limit: int | None = None,
    skip_on_error: bool = True,
) -> dict[str, int]:
    """Generate Run 4B mixed pseudo-label JSONL rows."""
    source = Path(input_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    table = load_codon_usage_table()
    counts = {
        "read": 0,
        "written": 0,
        "skipped_long": 0,
        "skipped_error": 0,
    }

    with output.open("w", encoding="utf-8") as handle:
        for protein_id, sequence in iter_protein_records(source, input_kind):
            if limit is not None and counts["written"] >= limit:
                break
            counts["read"] += 1
            protein = _normalize_protein(sequence)
            if len(protein) > max_length:
                protein = protein[:max_length]
                counts["skipped_long"] += 1
            try:
                row = build_mixed_teacher_row(
                    protein_id,
                    protein,
                    codon_weights=table.codon_weights,
                )
                handle.write(json.dumps(row, sort_keys=True) + "\n")
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
    rows = [
        line
        for line in Path(input_path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
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
    parser = argparse.ArgumentParser(description="Generate Run 4B mixed pseudo-label JSONL.")
    parser.add_argument("--input", required=True, help="Input protein/CDS FASTA or JSONL.")
    parser.add_argument("--output", default="data/training/run4b_mixed_pseudolabels.jsonl")
    parser.add_argument("--train-output", default=None)
    parser.add_argument("--eval-output", default=None)
    parser.add_argument("--train-split", type=float, default=0.9)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--max-length", type=int, default=512)
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument(
        "--input-kind",
        choices=["protein", "cds", "jsonl"],
        default="protein",
    )
    args = parser.parse_args()

    counts = generate_run4b_mixed_pseudolabels(
        input_path=args.input,
        output_path=args.output,
        input_kind=args.input_kind,
        max_length=args.max_length,
        limit=args.limit,
        skip_on_error=not args.fail_fast,
    )
    payload: dict[str, Any] = {
        "generation": counts,
        "available_v2_profiles": sorted(AVAILABLE_V2_PROFILES),
        "missing_v2_profiles": [
            profile
            for profile in ["local_gc_balanced", "motif_clean"]
            if profile not in AVAILABLE_V2_PROFILES
        ],
    }
    if args.train_output or args.eval_output:
        if not args.train_output or not args.eval_output:
            parser.error("--train-output and --eval-output must be provided together")
        payload["split"] = split_jsonl(
            args.output,
            args.train_output,
            args.eval_output,
            train_split=args.train_split,
            seed=args.seed,
        )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
