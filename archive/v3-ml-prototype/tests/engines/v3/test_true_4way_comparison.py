"""Tests for true 4-way sequence comparison."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EVALUATION = ROOT / "scripts" / "4_evaluation"
if str(EVALUATION) not in sys.path:
    sys.path.insert(0, str(EVALUATION))

from true_4way_comparison import (  # noqa: E402
    codon_identity,
    generate_feasibility_best,
    run_true_4way,
)
from factorforge.engines.v3.metrics import load_codon_usage_table  # noqa: E402


def _comparison_row(protein_id: str = "p1") -> dict[str, str]:
    return {
        "protein_id": protein_id,
        "v2_cai": "1.0",
        "v3_cai": "1.0",
        "v3_cai_delta_pct": "0.0000",
        "v2_gc_global": "50.0",
        "v3_gc_global": "50.0",
        "v2_local_gc_min": "50.0",
        "v2_local_gc_max": "50.0",
        "v3_local_gc_min": "50.0",
        "v3_local_gc_max": "50.0",
        "v3_amino_acid_identity": "1.0",
        "v3_validator_pass": "True",
        "v3_internal_stop_count": "0",
        "v3_invalid_codon_count": "0",
        "v3_forbidden_motif_count": "0",
        "v3_homopolymer_count": "0",
        "v3_repeat_count": "0",
    }


def _eval_row(protein_id: str = "p1") -> dict[str, object]:
    return {
        "protein_id": protein_id,
        "amino_acid_sequence": "MA",
        "protein_sequence": "MA",
        "sequence": "MA",
        "selected_teacher_type": "gc_target",
    }


def test_codon_identity_known_sequences() -> None:
    assert codon_identity("ATGGCCAAA", "ATGGCCAAA") == 1.0
    assert codon_identity("ATGGCCAAA", "ATGGCTAAA") == 2 / 3
    assert codon_identity("ATGGCC", "ATGGCCAAA") == 0.0


def test_feasibility_best_passes_validator() -> None:
    table = load_codon_usage_table()
    row = generate_feasibility_best(
        "p1",
        "MA",
        table.codon_weights,
        selected_teacher_type="gc_target",
    )

    assert row["candidate_type"] == "feasibility_best"
    assert row["feasibility_source"] == "dp_optimal_gc_40_55"
    assert row["validator_pass"] is True
    assert math.isfinite(float(row["multi_objective_score"]))


def test_run_true_4way_generates_all_candidate_types(tmp_path: Path) -> None:
    comparison = tmp_path / "comparison.csv"
    eval_file = tmp_path / "eval.jsonl"
    output_csv = tmp_path / "true_4way_comparison.csv"
    report = tmp_path / "true_4way_report.md"

    with comparison.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(_comparison_row().keys()))
        writer.writeheader()
        writer.writerow(_comparison_row())
    eval_file.write_text(json.dumps(_eval_row(), sort_keys=True) + "\n", encoding="utf-8")

    result = run_true_4way(
        comparison,
        eval_file,
        output_csv,
        report,
        progress_every=0,
    )

    assert result["rows"] == 4
    with output_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["candidate_type"] for row in rows} == {
        "v2_high_cai",
        "v2_gc_target",
        "v3_run4b_decoded",
        "feasibility_best",
    }
    assert [row for row in rows if row["candidate_type"] == "v3_run4b_decoded"][0][
        "candidate_status"
    ] == "redecode_required_no_local_sequence"
    assert "True 4-way Sequence Comparison" in report.read_text(encoding="utf-8")
