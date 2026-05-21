"""True 4-way sequence-level comparison for v3-alpha Run 4B outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2  # noqa: E402
from factorforge.engines.v3.metrics import load_codon_usage_table  # noqa: E402
from factorforge.ml.feasibility import analyze_feasibility  # noqa: E402
from factorforge.ml.metrics import (  # noqa: E402
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    detect_forbidden_motifs,
    detect_homopolymers,
    detect_repeats,
)
from factorforge.utils.validation import validate_candidate_sequence  # noqa: E402


GC_LOW = 40.0
GC_HIGH = 55.0
SCORE_WEIGHTS = {
    "cai_weight": 1.0,
    "gc_bonus": 0.1,
    "repeat_penalty": 0.05,
    "homopolymer_penalty": 0.05,
    "motif_penalty": 0.1,
}
V3_SEQUENCE_COLUMNS = (
    "v3_dna_sequence",
    "v3_sequence",
    "v3_decoded_sequence",
    "v3_decoded_cds",
    "v3_cds",
)
FIELDNAMES = [
    "protein_id",
    "candidate_type",
    "candidate_status",
    "dna_sequence",
    "sequence_available",
    "protein_length",
    "cai",
    "gc_global",
    "gc_in_40_55",
    "local_gc_min",
    "local_gc_max",
    "repeat_count",
    "homopolymer_count",
    "forbidden_motif_count",
    "multi_objective_score",
    "validator_pass",
    "amino_acid_identity",
    "internal_stop_count",
    "invalid_codon_count",
    "codon_identity_vs_v2_gc_target",
    "feasibility_source",
    "selected_teacher_type",
]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_jsonl_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def normalize_protein(row: dict[str, Any]) -> str:
    sequence = row.get("amino_acid_sequence") or row.get("protein_sequence") or row.get("sequence")
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("eval row is missing protein sequence metadata")
    return "".join(sequence.upper().split()).rstrip("*")


def multi_objective_score(
    *,
    cai: float,
    gc_pct: float,
    repeat_count: int,
    homopolymer_count: int,
    forbidden_motif_count: int,
) -> float:
    score = (
        SCORE_WEIGHTS["cai_weight"] * cai
        + SCORE_WEIGHTS["gc_bonus"] * (1.0 if GC_LOW <= gc_pct <= GC_HIGH else 0.0)
        - SCORE_WEIGHTS["repeat_penalty"] * repeat_count
        - SCORE_WEIGHTS["homopolymer_penalty"] * homopolymer_count
        - SCORE_WEIGHTS["motif_penalty"] * forbidden_motif_count
    )
    if not math.isfinite(score):
        raise ValueError("multi-objective score must be finite")
    return score


def codon_identity(seq_a: str, seq_b: str) -> float:
    """Return same-codon fraction over aligned complete codon positions."""
    codons_a = [seq_a[index : index + 3] for index in range(0, len(seq_a) - len(seq_a) % 3, 3)]
    codons_b = [seq_b[index : index + 3] for index in range(0, len(seq_b) - len(seq_b) % 3, 3)]
    if not codons_a or len(codons_a) != len(codons_b):
        return 0.0
    matches = sum(1 for left, right in zip(codons_a, codons_b) if left == right)
    return matches / len(codons_a)


def sequence_metric_row(
    *,
    protein_id: str,
    candidate_type: str,
    protein: str,
    dna_sequence: str,
    codon_weights: dict[str, float],
    selected_teacher_type: str,
    candidate_status: str = "sequence_available",
    feasibility_source: str = "",
    codon_identity_vs_v2_gc_target: float | None = None,
) -> dict[str, Any]:
    validator = validate_candidate_sequence(protein, dna_sequence)
    windows = calculate_gc_windows(dna_sequence)
    window_values = [float(window["gc"]) for window in windows]
    gc = calculate_gc(dna_sequence)
    repeat_count = len(detect_repeats(dna_sequence))
    homopolymer_count = len(detect_homopolymers(dna_sequence))
    forbidden_motif_count = len(detect_forbidden_motifs(dna_sequence, []))
    cai = calculate_cai(dna_sequence, codon_weights)
    return {
        "protein_id": protein_id,
        "candidate_type": candidate_type,
        "candidate_status": candidate_status,
        "dna_sequence": dna_sequence,
        "sequence_available": True,
        "protein_length": len(protein),
        "cai": cai,
        "gc_global": gc,
        "gc_in_40_55": GC_LOW <= gc <= GC_HIGH,
        "local_gc_min": min(window_values) if window_values else 0.0,
        "local_gc_max": max(window_values) if window_values else 0.0,
        "repeat_count": repeat_count,
        "homopolymer_count": homopolymer_count,
        "forbidden_motif_count": forbidden_motif_count,
        "multi_objective_score": multi_objective_score(
            cai=cai,
            gc_pct=gc,
            repeat_count=repeat_count,
            homopolymer_count=homopolymer_count,
            forbidden_motif_count=forbidden_motif_count,
        ),
        "validator_pass": validator["passed"],
        "amino_acid_identity": validator["amino_acid_identity"],
        "internal_stop_count": validator["internal_stop_count"],
        "invalid_codon_count": validator["invalid_codon_count"],
        "codon_identity_vs_v2_gc_target": (
            "" if codon_identity_vs_v2_gc_target is None else codon_identity_vs_v2_gc_target
        ),
        "feasibility_source": feasibility_source,
        "selected_teacher_type": selected_teacher_type,
    }


def v3_placeholder_row(
    *,
    comparison: dict[str, str],
    protein: str,
    selected_teacher_type: str,
    status: str,
) -> dict[str, Any]:
    gc = _float(comparison.get("v3_gc_global"))
    repeat_count = int(_float(comparison.get("v3_repeat_count")))
    homopolymer_count = int(_float(comparison.get("v3_homopolymer_count")))
    forbidden_motif_count = int(_float(comparison.get("v3_forbidden_motif_count")))
    cai = _float(comparison.get("v3_cai"))
    return {
        "protein_id": comparison["protein_id"],
        "candidate_type": "v3_run4b_decoded",
        "candidate_status": status,
        "dna_sequence": "",
        "sequence_available": False,
        "protein_length": len(protein),
        "cai": cai,
        "gc_global": gc,
        "gc_in_40_55": GC_LOW <= gc <= GC_HIGH,
        "local_gc_min": _float(comparison.get("v3_local_gc_min")),
        "local_gc_max": _float(comparison.get("v3_local_gc_max")),
        "repeat_count": repeat_count,
        "homopolymer_count": homopolymer_count,
        "forbidden_motif_count": forbidden_motif_count,
        "multi_objective_score": multi_objective_score(
            cai=cai,
            gc_pct=gc,
            repeat_count=repeat_count,
            homopolymer_count=homopolymer_count,
            forbidden_motif_count=forbidden_motif_count,
        ),
        "validator_pass": _bool(comparison.get("v3_validator_pass")),
        "amino_acid_identity": _float(comparison.get("v3_amino_acid_identity")),
        "internal_stop_count": int(_float(comparison.get("v3_internal_stop_count"))),
        "invalid_codon_count": int(_float(comparison.get("v3_invalid_codon_count"))),
        "codon_identity_vs_v2_gc_target": "",
        "feasibility_source": "",
        "selected_teacher_type": selected_teacher_type,
    }


def v3_sequence_from_comparison(row: dict[str, str]) -> str | None:
    for column in V3_SEQUENCE_COLUMNS:
        value = row.get(column)
        if value:
            return "".join(value.upper().replace("U", "T").split())
    return None


def generate_v2_candidate(
    protein_id: str,
    protein: str,
    profile: str,
    codon_weights: dict[str, float],
    selected_teacher_type: str,
) -> dict[str, Any]:
    result = optimize_with_v2(protein, options={"profile": profile, "scan_mode": "fast"})
    return sequence_metric_row(
        protein_id=protein_id,
        candidate_type=f"v2_{profile}",
        protein=protein,
        dna_sequence=result["dna_sequence"],
        codon_weights=codon_weights,
        selected_teacher_type=selected_teacher_type,
    )


def generate_feasibility_best(
    protein_id: str,
    protein: str,
    codon_weights: dict[str, float],
    selected_teacher_type: str,
) -> dict[str, Any]:
    feasibility = analyze_feasibility(
        protein,
        codon_weights,
        target_gc_low=GC_LOW,
        target_gc_high=GC_HIGH,
        gc_ranges=[(GC_LOW, GC_HIGH)],
    )
    range_result = feasibility["ranges"][f"{GC_LOW:g}-{GC_HIGH:g}"]
    candidate = range_result.get("best_candidate") if range_result.get("feasible") else None
    if not candidate or not isinstance(candidate.get("dna_sequence"), str):
        return {
            "protein_id": protein_id,
            "candidate_type": "feasibility_best",
            "candidate_status": "no_dp_candidate_under_gc_40_55",
            "dna_sequence": "",
            "sequence_available": False,
            "protein_length": len(protein),
            "cai": 0.0,
            "gc_global": 0.0,
            "gc_in_40_55": False,
            "local_gc_min": 0.0,
            "local_gc_max": 0.0,
            "repeat_count": 0,
            "homopolymer_count": 0,
            "forbidden_motif_count": 0,
            "multi_objective_score": 0.0,
            "validator_pass": False,
            "amino_acid_identity": 0.0,
            "internal_stop_count": 0,
            "invalid_codon_count": 0,
            "codon_identity_vs_v2_gc_target": "",
            "feasibility_source": "dp_no_candidate",
            "selected_teacher_type": selected_teacher_type,
        }
    return sequence_metric_row(
        protein_id=protein_id,
        candidate_type="feasibility_best",
        protein=protein,
        dna_sequence=candidate["dna_sequence"],
        codon_weights=codon_weights,
        selected_teacher_type=selected_teacher_type,
        candidate_status="sequence_available",
        feasibility_source="dp_optimal_gc_40_55",
    )


def build_true_4way_rows(
    comparison_rows: list[dict[str, str]],
    eval_rows: list[dict[str, Any]],
    codon_weights: dict[str, float] | None = None,
    progress_every: int = 250,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    weights = codon_weights or load_codon_usage_table().codon_weights
    eval_by_id = {str(row["protein_id"]): row for row in eval_rows}
    present_columns = set(comparison_rows[0]) if comparison_rows else set()
    v3_sequence_available = any(column in present_columns for column in V3_SEQUENCE_COLUMNS)
    v3_status = "sequence_available" if v3_sequence_available else "redecode_required_no_local_sequence"

    rows: list[dict[str, Any]] = []
    identity_values: list[float] = []
    start = time.time()
    for index, comparison in enumerate(comparison_rows, start=1):
        protein_id = comparison["protein_id"]
        eval_row = eval_by_id[protein_id]
        protein = normalize_protein(eval_row)
        selected_teacher_type = str(eval_row.get("selected_teacher_type", ""))

        high_cai = generate_v2_candidate(
            protein_id, protein, "high_cai", weights, selected_teacher_type
        )
        gc_target = generate_v2_candidate(
            protein_id, protein, "gc_target", weights, selected_teacher_type
        )
        v3_sequence = v3_sequence_from_comparison(comparison)
        if v3_sequence:
            identity = codon_identity(v3_sequence, gc_target["dna_sequence"])
            v3 = sequence_metric_row(
                protein_id=protein_id,
                candidate_type="v3_run4b_decoded",
                protein=protein,
                dna_sequence=v3_sequence,
                codon_weights=weights,
                selected_teacher_type=selected_teacher_type,
                candidate_status="sequence_available",
                codon_identity_vs_v2_gc_target=identity,
            )
            identity_values.append(identity)
        else:
            v3 = v3_placeholder_row(
                comparison=comparison,
                protein=protein,
                selected_teacher_type=selected_teacher_type,
                status=v3_status,
            )
        feasibility_best = generate_feasibility_best(
            protein_id, protein, weights, selected_teacher_type
        )
        rows.extend([high_cai, gc_target, v3, feasibility_best])

        if progress_every > 0 and (
            index == 1 or index % progress_every == 0 or index == len(comparison_rows)
        ):
            elapsed = time.time() - start
            rate = index / elapsed if elapsed > 0 else 0.0
            print(
                f"true_4way {index}/{len(comparison_rows)} proteins | "
                f"{rate:.2f} proteins/sec",
                flush=True,
            )

    diagnostics = {
        "present_columns": sorted(present_columns),
        "missing_v3_sequence_columns": [
            column for column in V3_SEQUENCE_COLUMNS if column not in present_columns
        ],
        "v3_sequence_status": v3_status,
        "v3_vs_gc_target_identity_mean": (
            statistics.fmean(identity_values) if identity_values else None
        ),
    }
    return rows, diagnostics


def summary_stats(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["candidate_type"])].append(row)
    summary: dict[str, dict[str, float]] = {}
    for candidate_type, items in grouped.items():
        summary[candidate_type] = {
            "n": float(len(items)),
            "sequence_available_rate": statistics.fmean(
                1.0 if _bool(row["sequence_available"]) else 0.0 for row in items
            ),
            "cai_mean": statistics.fmean(_float(row["cai"]) for row in items),
            "gc_mean": statistics.fmean(_float(row["gc_global"]) for row in items),
            "repeat_mean": statistics.fmean(_float(row["repeat_count"]) for row in items),
            "homopolymer_mean": statistics.fmean(
                _float(row["homopolymer_count"]) for row in items
            ),
            "motif_mean": statistics.fmean(_float(row["forbidden_motif_count"]) for row in items),
            "score_mean": statistics.fmean(
                _float(row["multi_objective_score"]) for row in items
            ),
            "validator_pass_rate": statistics.fmean(
                1.0 if _bool(row["validator_pass"]) else 0.0 for row in items
            ),
        }
    return summary


def _fmt(value: float | None) -> str:
    return "N/A" if value is None else f"{value:.4f}"


def summary_table(summary: dict[str, dict[str, float]]) -> str:
    lines = [
        "| candidate_type | n | sequence_available | CAI mean | GC mean | repeat mean | homopolymer mean | motif mean | score mean | validator pass |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for candidate_type in [
        "v2_high_cai",
        "v2_gc_target",
        "v3_run4b_decoded",
        "feasibility_best",
    ]:
        stats = summary.get(candidate_type)
        if not stats:
            continue
        lines.append(
            "| "
            + " | ".join(
                [
                    candidate_type,
                    str(int(stats["n"])),
                    _fmt(stats["sequence_available_rate"]),
                    _fmt(stats["cai_mean"]),
                    _fmt(stats["gc_mean"]),
                    _fmt(stats["repeat_mean"]),
                    _fmt(stats["homopolymer_mean"]),
                    _fmt(stats["motif_mean"]),
                    _fmt(stats["score_mean"]),
                    _fmt(stats["validator_pass_rate"]),
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def build_report(
    rows: list[dict[str, Any]],
    diagnostics: dict[str, Any],
    comparison_path: Path,
    eval_path: Path,
    output_csv: Path,
) -> str:
    summary = summary_stats(rows)
    status_counts = Counter(str(row["candidate_status"]) for row in rows)
    feasibility_counts = Counter(
        str(row["feasibility_source"])
        for row in rows
        if row["candidate_type"] == "feasibility_best"
    )
    identity_mean = diagnostics["v3_vs_gc_target_identity_mean"]
    v3_status = diagnostics["v3_sequence_status"]
    if identity_mean is None:
        independence = (
            "Cannot determine from local artifacts because v3 decoded CDS sequences are not "
            "stored in alpha_run2 outputs. Re-decode with the Run4B checkpoint is required."
        )
    elif identity_mean >= 0.99:
        independence = "v3 is essentially a codon-level copy of v2 gc_target on this eval set."
    else:
        independence = "v3 has a distinct codon-selection pattern from v2 gc_target on this eval set."

    return "\n".join(
        [
            "# True 4-way Sequence Comparison",
            "",
            "## Inputs",
            "",
            f"- comparison_csv: `{comparison_path}`",
            f"- eval_jsonl: `{eval_path}`",
            f"- output_csv: `{output_csv}`",
            f"- proteins_evaluated: {int(summary.get('v2_high_cai', {}).get('n', 0.0))}",
            "",
            "## Source Column Check",
            "",
            f"- present_columns: {', '.join(diagnostics['present_columns'])}",
            f"- missing_v3_sequence_columns: {', '.join(diagnostics['missing_v3_sequence_columns'])}",
            f"- v3_sequence_status: {v3_status}",
            "",
            "Candidate C uses stored v3 decoded CDS sequences only if such a column exists. "
            "The current local alpha_run2 comparison stores v3 metrics but not decoded CDS, "
            "so v3 rows are metric placeholders and require re-decode to support codon identity.",
            "",
            "## 4-way Comparison Summary",
            "",
            summary_table(summary),
            "",
            "## v3 vs v2 gc_target Codon Identity",
            "",
            f"- mean codon identity: {_fmt(identity_mean)}",
            f"- interpretation: {independence}",
            "",
            "## Feasibility Candidate",
            "",
            f"- feasibility_source_counts: {dict(sorted(feasibility_counts.items()))}",
            "- `dp_optimal_gc_40_55` means `analyze_feasibility()` selected the maximum-CAI sequence under global GC 40-55%.",
            "",
            "## Candidate Status Counts",
            "",
            "```json",
            json.dumps(dict(sorted(status_counts.items())), indent=2, sort_keys=True),
            "```",
            "",
            "## Honest Conclusion",
            "",
            independence,
            "",
            "This is an in-silico sequence comparison. It does not change v2 production behavior "
            "and does not support yield or expression claims.",
            "",
        ]
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def run_true_4way(
    comparison_path: Path,
    eval_path: Path,
    output_csv: Path,
    report_path: Path,
    limit: int | None = None,
    progress_every: int = 250,
) -> dict[str, Any]:
    comparison_rows = read_csv_rows(comparison_path)
    eval_rows = read_jsonl_rows(eval_path)
    if limit is not None:
        comparison_rows = comparison_rows[:limit]
    rows, diagnostics = build_true_4way_rows(
        comparison_rows,
        eval_rows,
        progress_every=progress_every,
    )
    write_csv(output_csv, rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        build_report(rows, diagnostics, comparison_path, eval_path, output_csv),
        encoding="utf-8",
    )
    feasibility_counts = Counter(
        str(row["feasibility_source"])
        for row in rows
        if row["candidate_type"] == "feasibility_best"
    )
    return {
        "proteins": len(comparison_rows),
        "rows": len(rows),
        "output_csv": str(output_csv),
        "report": str(report_path),
        "v3_sequence_status": diagnostics["v3_sequence_status"],
        "v3_vs_gc_target_identity_mean": diagnostics["v3_vs_gc_target_identity_mean"],
        "feasibility_source_counts": dict(sorted(feasibility_counts.items())),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run true 4-way v3-alpha comparison.")
    parser.add_argument(
        "--comparison",
        type=Path,
        default=ROOT / "experiments" / "results" / "alpha_run2" / "alpha_run2_comparison.csv",
    )
    parser.add_argument(
        "--eval-file",
        type=Path,
        default=ROOT / "data" / "training" / "run4b_pseudolabels_eval.jsonl",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=ROOT / "experiments" / "results" / "alpha_run2" / "true_4way_comparison.csv",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "experiments" / "results" / "alpha_run2" / "true_4way_report.md",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=250)
    args = parser.parse_args()
    print(
        json.dumps(
            run_true_4way(
                comparison_path=args.comparison,
                eval_path=args.eval_file,
                output_csv=args.output_csv,
                report_path=args.report,
                limit=args.limit,
                progress_every=args.progress_every,
            ),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
