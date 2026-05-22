"""Build v3-alpha evidence/reranking evaluation tables from Run 4B outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
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


DEFAULT_WEIGHTS = {
    "cai_weight": 1.0,
    "gc_bonus": 0.1,
    "repeat_penalty": 0.05,
    "homopolymer_penalty": 0.05,
    "motif_penalty": 0.1,
}
GC_LOW = 40.0
GC_HIGH = 55.0
REQUIRED_EXTENDED_COLUMNS = {
    "protein_id",
    "candidate_type",
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
}


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


def multi_objective_score(
    *,
    cai: float,
    gc_pct: float,
    repeat_count: int,
    homopolymer_count: int,
    forbidden_motif_count: int,
    weights: dict[str, float] | None = None,
) -> float:
    """Score one candidate using the v3-alpha evidence weights."""
    w = {**DEFAULT_WEIGHTS, **(weights or {})}
    score = (
        w["cai_weight"] * cai
        + w["gc_bonus"] * (1.0 if GC_LOW <= gc_pct <= GC_HIGH else 0.0)
        - w["repeat_penalty"] * repeat_count
        - w["homopolymer_penalty"] * homopolymer_count
        - w["motif_penalty"] * forbidden_motif_count
    )
    if not math.isfinite(score):
        raise ValueError("multi-objective score must be finite")
    return score


def _protein_sequence(row: dict[str, Any]) -> str:
    sequence = row.get("amino_acid_sequence") or row.get("protein_sequence") or row.get("sequence")
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("row is missing amino acid sequence metadata")
    return "".join(sequence.upper().split()).rstrip("*")


def _sequence_metric_row(
    *,
    protein_id: str,
    candidate_type: str,
    protein: str,
    dna_sequence: str,
    codon_weights: dict[str, float],
    selected_teacher_type: str = "",
    source_sequence_available: bool = True,
) -> dict[str, Any]:
    validator = validate_candidate_sequence(protein, dna_sequence)
    windows = calculate_gc_windows(dna_sequence)
    window_values = [float(window["gc"]) for window in windows]
    gc = calculate_gc(dna_sequence)
    repeat_count = len(detect_repeats(dna_sequence))
    homopolymer_count = len(detect_homopolymers(dna_sequence))
    forbidden_motif_count = len(detect_forbidden_motifs(dna_sequence, []))
    cai = calculate_cai(dna_sequence, codon_weights)
    score = multi_objective_score(
        cai=cai,
        gc_pct=gc,
        repeat_count=repeat_count,
        homopolymer_count=homopolymer_count,
        forbidden_motif_count=forbidden_motif_count,
    )
    return {
        "protein_id": protein_id,
        "candidate_type": candidate_type,
        "selected_teacher_type": selected_teacher_type,
        "protein_length": len(protein),
        "candidate_available": True,
        "source_sequence_available": source_sequence_available,
        "feasibility_source": "",
        "cai": cai,
        "gc_global": gc,
        "gc_in_40_55": GC_LOW <= gc <= GC_HIGH,
        "local_gc_min": min(window_values) if window_values else 0.0,
        "local_gc_max": max(window_values) if window_values else 0.0,
        "repeat_count": repeat_count,
        "homopolymer_count": homopolymer_count,
        "forbidden_motif_count": forbidden_motif_count,
        "multi_objective_score": score,
        "amino_acid_identity": validator["amino_acid_identity"],
        "validator_pass": validator["passed"],
        "internal_stop_count": validator["internal_stop_count"],
        "invalid_codon_count": validator["invalid_codon_count"],
        "gc_window_outlier_count": validator["gc_window_outlier_count"],
    }


def _v2_candidate_row(
    protein_id: str,
    candidate_type: str,
    profile: str,
    protein: str,
    codon_weights: dict[str, float],
    selected_teacher_type: str,
) -> dict[str, Any]:
    result = optimize_with_v2(protein, options={"profile": profile, "scan_mode": "fast"})
    return _sequence_metric_row(
        protein_id=protein_id,
        candidate_type=candidate_type,
        protein=protein,
        dna_sequence=result["dna_sequence"],
        codon_weights=codon_weights,
        selected_teacher_type=selected_teacher_type,
    )


def _metric_summary_row(
    *,
    protein_id: str,
    candidate_type: str,
    protein: str,
    summary: dict[str, Any],
    selected_teacher_type: str,
    feasibility_source: str = "",
) -> dict[str, Any]:
    gc = _float(summary.get("gc"))
    cai = _float(summary.get("cai"))
    repeat_count = int(_float(summary.get("repeat_count")))
    homopolymer_count = int(_float(summary.get("homopolymer_count")))
    forbidden_motif_count = int(_float(summary.get("forbidden_motif_count")))
    validator = summary.get("validator") if isinstance(summary.get("validator"), dict) else {}
    score = multi_objective_score(
        cai=cai,
        gc_pct=gc,
        repeat_count=repeat_count,
        homopolymer_count=homopolymer_count,
        forbidden_motif_count=forbidden_motif_count,
    )
    return {
        "protein_id": protein_id,
        "candidate_type": candidate_type,
        "selected_teacher_type": selected_teacher_type,
        "protein_length": len(protein),
        "candidate_available": True,
        "source_sequence_available": False,
        "feasibility_source": feasibility_source,
        "cai": cai,
        "gc_global": gc,
        "gc_in_40_55": bool(summary.get("gc_in_40_55", GC_LOW <= gc <= GC_HIGH)),
        "local_gc_min": _float(summary.get("gc_window_min")),
        "local_gc_max": _float(summary.get("gc_window_max")),
        "repeat_count": repeat_count,
        "homopolymer_count": homopolymer_count,
        "forbidden_motif_count": forbidden_motif_count,
        "multi_objective_score": score,
        "amino_acid_identity": _float(validator.get("amino_acid_identity"), 1.0),
        "validator_pass": bool(validator.get("passed", True)),
        "internal_stop_count": int(_float(validator.get("internal_stop_count"))),
        "invalid_codon_count": int(_float(validator.get("invalid_codon_count"))),
        "gc_window_outlier_count": int(_float(validator.get("gc_window_outlier_count"))),
    }


def _candidate_metric_from_eval(
    eval_row: dict[str, Any],
    candidate_name: str,
) -> dict[str, Any] | None:
    candidate_metrics = eval_row.get("candidate_metrics")
    if not isinstance(candidate_metrics, dict):
        return None
    candidate = candidate_metrics.get(candidate_name)
    return candidate if isinstance(candidate, dict) else None


def _feasibility_candidate_row(
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
    range_result = feasibility["ranges"].get(f"{GC_LOW:g}-{GC_HIGH:g}", {})
    candidate = (range_result.get("best_candidate") or {}) if range_result else {}
    dna_sequence = candidate.get("dna_sequence")
    if not isinstance(dna_sequence, str) or not dna_sequence:
        return {
            "protein_id": protein_id,
            "candidate_type": "feasibility_best",
            "selected_teacher_type": selected_teacher_type,
            "protein_length": len(protein),
            "candidate_available": False,
            "source_sequence_available": False,
            "feasibility_source": "unavailable",
            "cai": 0.0,
            "gc_global": 0.0,
            "gc_in_40_55": False,
            "local_gc_min": 0.0,
            "local_gc_max": 0.0,
            "repeat_count": 0,
            "homopolymer_count": 0,
            "forbidden_motif_count": 0,
            "multi_objective_score": 0.0,
            "amino_acid_identity": 0.0,
            "validator_pass": False,
            "internal_stop_count": 0,
            "invalid_codon_count": 0,
            "gc_window_outlier_count": 0,
        }
    return _sequence_metric_row(
        protein_id=protein_id,
        candidate_type="feasibility_best",
        protein=protein,
        dna_sequence=dna_sequence,
        codon_weights=codon_weights,
        selected_teacher_type=selected_teacher_type,
    )


def _v3_decoded_metric_row(
    comparison: dict[str, str],
    eval_row: dict[str, Any],
) -> dict[str, Any]:
    gc = _float(comparison.get("v3_gc_global"))
    repeat_count = int(_float(comparison.get("v3_repeat_count")))
    homopolymer_count = int(_float(comparison.get("v3_homopolymer_count")))
    forbidden_motif_count = int(_float(comparison.get("v3_forbidden_motif_count")))
    cai = _float(comparison.get("v3_cai"))
    score = multi_objective_score(
        cai=cai,
        gc_pct=gc,
        repeat_count=repeat_count,
        homopolymer_count=homopolymer_count,
        forbidden_motif_count=forbidden_motif_count,
    )
    protein = _protein_sequence(eval_row)
    return {
        "protein_id": comparison["protein_id"],
        "candidate_type": "v3_run4b_decoded",
        "selected_teacher_type": str(eval_row.get("selected_teacher_type", "")),
        "protein_length": len(protein),
        "candidate_available": True,
        "source_sequence_available": False,
        "feasibility_source": "",
        "cai": cai,
        "gc_global": gc,
        "gc_in_40_55": GC_LOW <= gc <= GC_HIGH,
        "local_gc_min": _float(comparison.get("v3_local_gc_min")),
        "local_gc_max": _float(comparison.get("v3_local_gc_max")),
        "repeat_count": repeat_count,
        "homopolymer_count": homopolymer_count,
        "forbidden_motif_count": forbidden_motif_count,
        "multi_objective_score": score,
        "amino_acid_identity": _float(comparison.get("v3_amino_acid_identity")),
        "validator_pass": _bool(comparison.get("v3_validator_pass")),
        "internal_stop_count": int(_float(comparison.get("v3_internal_stop_count"))),
        "invalid_codon_count": int(_float(comparison.get("v3_invalid_codon_count"))),
        "gc_window_outlier_count": "",
    }


def build_extended_rows(
    comparison_rows: list[dict[str, str]],
    eval_rows: list[dict[str, Any]],
    codon_weights: dict[str, float] | None = None,
    compute_missing_feasibility: bool = False,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build long-form 4-way candidate rows and source-column diagnostics."""
    weights = codon_weights or load_codon_usage_table().codon_weights
    eval_by_id = {str(row["protein_id"]): row for row in eval_rows}
    present_columns = set(comparison_rows[0].keys()) if comparison_rows else set()
    expected_columns = {
        "v2_repeat_count",
        "v2_homopolymer_count",
        "v2_forbidden_motif_count",
        "v2_local_gc_min",
        "v2_local_gc_max",
        "v3_repeat_count",
        "v3_homopolymer_count",
        "v3_forbidden_motif_count",
        "v3_local_gc_min",
        "v3_local_gc_max",
        "v2_sequence",
        "v3_sequence",
    }
    diagnostics = {
        "present_columns": sorted(present_columns),
        "missing_columns": sorted(expected_columns - present_columns),
    }

    rows: list[dict[str, Any]] = []
    for comparison in comparison_rows:
        protein_id = comparison["protein_id"]
        eval_row = eval_by_id.get(protein_id)
        if eval_row is None:
            raise KeyError(f"Missing eval row for protein_id={protein_id}")
        protein = _protein_sequence(eval_row)
        selected_teacher_type = str(eval_row.get("selected_teacher_type", ""))
        high_cai_summary = _candidate_metric_from_eval(eval_row, "high_cai")
        if high_cai_summary is None:
            rows.append(
                _v2_candidate_row(
                    protein_id,
                    "v2_high_cai",
                    "high_cai",
                    protein,
                    weights,
                    selected_teacher_type,
                )
            )
        else:
            rows.append(
                _metric_summary_row(
                    protein_id=protein_id,
                    candidate_type="v2_high_cai",
                    protein=protein,
                    summary=high_cai_summary,
                    selected_teacher_type=selected_teacher_type,
                )
            )

        gc_target_summary = _candidate_metric_from_eval(eval_row, "gc_target")
        if gc_target_summary is None:
            rows.append(
                _v2_candidate_row(
                    protein_id,
                    "v2_gc_target",
                    "gc_target",
                    protein,
                    weights,
                    selected_teacher_type,
                )
            )
        else:
            rows.append(
                _metric_summary_row(
                    protein_id=protein_id,
                    candidate_type="v2_gc_target",
                    protein=protein,
                    summary=gc_target_summary,
                    selected_teacher_type=selected_teacher_type,
                )
            )
        rows.append(_v3_decoded_metric_row(comparison, eval_row))
        feasibility_summary = _candidate_metric_from_eval(eval_row, "feasibility_40_55")
        if feasibility_summary is None and compute_missing_feasibility:
            rows.append(_feasibility_candidate_row(protein_id, protein, weights, selected_teacher_type))
        elif feasibility_summary is not None:
            rows.append(
                _metric_summary_row(
                    protein_id=protein_id,
                    candidate_type="feasibility_best",
                    protein=protein,
                    summary=feasibility_summary,
                    selected_teacher_type=selected_teacher_type,
                    feasibility_source="analyze_feasibility_jsonl",
                )
            )
        elif gc_target_summary is not None:
            rows.append(
                _metric_summary_row(
                    protein_id=protein_id,
                    candidate_type="feasibility_best",
                    protein=protein,
                    summary=gc_target_summary,
                    selected_teacher_type=selected_teacher_type,
                    feasibility_source="gc_target_fallback_exact_not_in_jsonl",
                )
            )
        else:
            rows.append(_feasibility_candidate_row(protein_id, protein, weights, selected_teacher_type))
    return rows, diagnostics


def _format_number(value: float) -> str:
    return f"{value:.4f}"


def _summary_stats(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["candidate_type"])].append(row)
    summary: dict[str, dict[str, float]] = {}
    for candidate_type, items in grouped.items():
        scores = [float(row["multi_objective_score"]) for row in items]
        cai = [float(row["cai"]) for row in items]
        gc = [float(row["gc_global"]) for row in items]
        summary[candidate_type] = {
            "n": float(len(items)),
            "score_mean": statistics.fmean(scores),
            "score_median": statistics.median(scores),
            "score_min": min(scores),
            "score_max": max(scores),
            "cai_mean": statistics.fmean(cai),
            "gc_mean": statistics.fmean(gc),
            "gc_in_range_rate": statistics.fmean(
                1.0 if bool(row["gc_in_40_55"]) else 0.0 for row in items
            ),
            "validator_pass_rate": statistics.fmean(
                1.0 if bool(row["validator_pass"]) else 0.0 for row in items
            ),
        }
    return summary


def _score_by_protein(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    scores: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        scores[str(row["protein_id"])][str(row["candidate_type"])] = float(
            row["multi_objective_score"]
        )
    return scores


def _v3_advantage(rows: list[dict[str, Any]]) -> dict[str, Any]:
    score_map = _score_by_protein(rows)
    v3_gt_high = [
        protein_id
        for protein_id, scores in score_map.items()
        if scores.get("v3_run4b_decoded", -math.inf) > scores.get("v2_high_cai", math.inf)
    ]
    v3_gt_gc = [
        protein_id
        for protein_id, scores in score_map.items()
        if scores.get("v3_run4b_decoded", -math.inf) > scores.get("v2_gc_target", math.inf)
    ]
    by_id = {str(row["protein_id"]): row for row in rows if row["candidate_type"] == "v3_run4b_decoded"}

    def group_characteristics(ids: list[str]) -> dict[str, Any]:
        items = [by_id[item] for item in ids if item in by_id]
        selected = Counter(str(row["selected_teacher_type"]) for row in items)
        lengths = [float(row["protein_length"]) for row in items]
        return {
            "count": len(items),
            "protein_length_mean": statistics.fmean(lengths) if lengths else 0.0,
            "selected_teacher_type_counts": dict(sorted(selected.items())),
        }

    return {
        "v3_gt_v2_high_cai_count": len(v3_gt_high),
        "v3_gt_v2_gc_target_count": len(v3_gt_gc),
        "v3_gt_v2_high_cai_ids": v3_gt_high[:20],
        "v3_gt_v2_gc_target_ids": v3_gt_gc[:20],
        "v3_gt_v2_high_cai_characteristics": group_characteristics(v3_gt_high),
        "v3_gt_v2_gc_target_characteristics": group_characteristics(v3_gt_gc),
    }


def _summary_table_markdown(summary: dict[str, dict[str, float]]) -> str:
    lines = [
        "| candidate_type | n | score_mean | score_median | score_min | score_max | cai_mean | gc_mean | gc_in_range_rate | validator_pass_rate |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for candidate_type in sorted(summary):
        stats = summary[candidate_type]
        lines.append(
            "| "
            + " | ".join(
                [
                    candidate_type,
                    str(int(stats["n"])),
                    _format_number(stats["score_mean"]),
                    _format_number(stats["score_median"]),
                    _format_number(stats["score_min"]),
                    _format_number(stats["score_max"]),
                    _format_number(stats["cai_mean"]),
                    _format_number(stats["gc_mean"]),
                    _format_number(stats["gc_in_range_rate"]),
                    _format_number(stats["validator_pass_rate"]),
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def _conclusion(summary: dict[str, dict[str, float]], advantage: dict[str, Any]) -> str:
    v3 = summary.get("v3_run4b_decoded", {})
    high = summary.get("v2_high_cai", {})
    gc = summary.get("v2_gc_target", {})
    if not v3 or not high or not gc:
        return "no advantage found"
    total = int(v3.get("n", 0))
    if total and (
        advantage["v3_gt_v2_high_cai_count"] > total * 0.5
        and advantage["v3_gt_v2_gc_target_count"] > total * 0.5
    ):
        return "v3 is better"
    mean_gap_high = float(v3["score_mean"]) - float(high["score_mean"])
    mean_gap_gc = float(v3["score_mean"]) - float(gc["score_mean"])
    if abs(mean_gap_high) <= 0.02 or abs(mean_gap_gc) <= 0.02:
        return "v3 is comparable"
    return "no advantage found"


def build_report(
    *,
    rows: list[dict[str, Any]],
    diagnostics: dict[str, Any],
    comparison_path: Path,
    eval_path: Path,
    output_csv: Path,
) -> str:
    summary = _summary_stats(rows)
    advantage = _v3_advantage(rows)
    conclusion = _conclusion(summary, advantage)
    total_proteins = int(summary.get("v3_run4b_decoded", {}).get("n", 0.0))
    feasibility_sources = Counter(
        str(row.get("feasibility_source", ""))
        for row in rows
        if row.get("candidate_type") == "feasibility_best"
    )

    return "\n".join(
        [
            "# v3-alpha Evidence/Reranking Evaluation",
            "",
            "## Inputs",
            "",
            f"- comparison_csv: `{comparison_path}`",
            f"- eval_jsonl: `{eval_path}`",
            f"- extended_csv: `{output_csv}`",
            f"- proteins_evaluated: {total_proteins}",
            "",
            "## Source CSV Columns",
            "",
            f"- present: {', '.join(diagnostics['present_columns'])}",
            f"- missing: {', '.join(diagnostics['missing_columns']) or 'none'}",
            "",
            "The Run 4B comparison CSV does not include decoded CDS sequence columns. "
            "The v3 Run4B row therefore uses the decoded metrics already recorded in "
            "`run4_comparison.csv`; v2 high_cai, v2 gc_target, and feasibility_best "
            "rows use the sequence-level metrics recorded during Run4B teacher generation.",
            "",
            "For feasibility_best, exact `analyze_feasibility` rows are used when they "
            "were recorded in the Run4B eval JSONL. Rows without an exact feasibility "
            "candidate in JSONL use the valid gc_target candidate as a fallback and are "
            "marked in `feasibility_source`.",
            "",
            f"- feasibility_source_counts: {dict(sorted(feasibility_sources.items()))}",
            "",
            "## 4-way Comparison Summary",
            "",
            _summary_table_markdown(summary),
            "",
            "## Multi-objective Score",
            "",
            "Score formula: `cai + 0.1 * GC_in_40_55 - 0.05 * repeats "
            "- 0.05 * homopolymers - 0.1 * forbidden_motifs`.",
            "",
            "## v3 Advantage Cases",
            "",
            f"- v3 score > v2 high_cai: {advantage['v3_gt_v2_high_cai_count']} / {total_proteins}",
            f"- v3 score > v2 gc_target: {advantage['v3_gt_v2_gc_target_count']} / {total_proteins}",
            f"- first v3 > v2 high_cai IDs: {advantage['v3_gt_v2_high_cai_ids']}",
            f"- first v3 > v2 gc_target IDs: {advantage['v3_gt_v2_gc_target_ids']}",
            "",
            "## Protein Characteristics Associated With v3 Advantage",
            "",
            "v3 > v2 high_cai:",
            "",
            "```json",
            json.dumps(
                advantage["v3_gt_v2_high_cai_characteristics"],
                indent=2,
                sort_keys=True,
            ),
            "```",
            "",
            "v3 > v2 gc_target:",
            "",
            "```json",
            json.dumps(
                advantage["v3_gt_v2_gc_target_characteristics"],
                indent=2,
                sort_keys=True,
            ),
            "```",
            "",
            "## Honest Conclusion",
            "",
            conclusion,
            "",
            "This is an in-silico candidate evidence comparison only. It does not change v2 "
            "production behavior and does not support yield or expression claims.",
            "",
        ]
    )


def write_extended_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No rows to write")
    missing = REQUIRED_EXTENDED_COLUMNS - set(rows[0].keys())
    if missing:
        raise ValueError(f"extended rows missing required columns: {sorted(missing)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_evaluation(
    comparison_path: Path,
    eval_path: Path,
    output_csv: Path,
    report_path: Path,
    limit: int | None = None,
    compute_missing_feasibility: bool = False,
) -> dict[str, Any]:
    comparison_rows = read_csv_rows(comparison_path)
    eval_rows = read_jsonl_rows(eval_path)
    if limit is not None:
        comparison_rows = comparison_rows[:limit]
    rows, diagnostics = build_extended_rows(
        comparison_rows,
        eval_rows,
        compute_missing_feasibility=compute_missing_feasibility,
    )
    write_extended_csv(output_csv, rows)
    report = build_report(
        rows=rows,
        diagnostics=diagnostics,
        comparison_path=comparison_path,
        eval_path=eval_path,
        output_csv=output_csv,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    summary = _summary_stats(rows)
    advantage = _v3_advantage(rows)
    return {
        "rows": len(rows),
        "proteins": len(comparison_rows),
        "extended_csv": str(output_csv),
        "report": str(report_path),
        "conclusion": _conclusion(summary, advantage),
        "v3_gt_v2_high_cai_count": advantage["v3_gt_v2_high_cai_count"],
        "v3_gt_v2_gc_target_count": advantage["v3_gt_v2_gc_target_count"],
    }


def _default_comparison_path() -> Path:
    return ROOT / "experiments" / "results" / "alpha_run2" / "alpha_run2_comparison.csv"


def main() -> int:
    parser = argparse.ArgumentParser(description="Run v3-alpha evidence evaluation.")
    parser.add_argument(
        "--comparison",
        type=Path,
        default=_default_comparison_path(),
    )
    parser.add_argument(
        "--eval-file",
        type=Path,
        default=ROOT / "data" / "training" / "run4b_pseudolabels_eval.jsonl",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=ROOT
        / "experiments"
        / "results"
        / "alpha_run2"
        / "alpha_run2_comparison_extended.csv",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT
        / "experiments"
        / "results"
        / "alpha_run2"
        / "v3alpha_evidence_report.md",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--compute-missing-feasibility",
        action="store_true",
        help="Run exact analyze_feasibility for rows not already carrying feasibility metrics.",
    )
    args = parser.parse_args()
    print(
        json.dumps(
            run_evaluation(
                comparison_path=args.comparison,
                eval_path=args.eval_file,
                output_csv=args.output_csv,
                report_path=args.report,
                limit=args.limit,
                compute_missing_feasibility=args.compute_missing_feasibility,
            ),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
