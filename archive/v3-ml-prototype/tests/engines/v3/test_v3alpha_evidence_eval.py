"""Tests for v3-alpha evidence evaluation outputs."""

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

from v3alpha_evidence_eval import (  # noqa: E402
    REQUIRED_EXTENDED_COLUMNS,
    build_extended_rows,
    multi_objective_score,
    run_evaluation,
)


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
        "dna_sequence": "ATGGCC",
        "selected_teacher_type": "gc_target",
    }


def test_multi_objective_score_is_finite() -> None:
    score = multi_objective_score(
        cai=0.8,
        gc_pct=42.0,
        repeat_count=1,
        homopolymer_count=0,
        forbidden_motif_count=0,
    )

    assert math.isfinite(score)
    assert score == 0.85


def test_extended_rows_have_required_columns_and_finite_scores() -> None:
    rows, diagnostics = build_extended_rows([_comparison_row()], [_eval_row()])

    assert len(rows) == 4
    assert REQUIRED_EXTENDED_COLUMNS <= set(rows[0])
    assert "v2_repeat_count" in diagnostics["missing_columns"]
    assert "v3_sequence" in diagnostics["missing_columns"]
    assert all(math.isfinite(float(row["multi_objective_score"])) for row in rows)


def test_v2_candidates_pass_validator() -> None:
    rows, _ = build_extended_rows([_comparison_row()], [_eval_row()])
    v2_rows = [row for row in rows if str(row["candidate_type"]).startswith("v2_")]

    assert len(v2_rows) == 2
    assert all(row["validator_pass"] is True for row in v2_rows)
    assert all(float(row["amino_acid_identity"]) == 1.0 for row in v2_rows)


def test_run_evaluation_writes_extended_csv_and_report(tmp_path: Path) -> None:
    comparison = tmp_path / "run4_comparison.csv"
    eval_file = tmp_path / "eval.jsonl"
    output_csv = tmp_path / "run4_comparison_extended.csv"
    report = tmp_path / "v3alpha_evidence_report.md"

    with comparison.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(_comparison_row().keys()))
        writer.writeheader()
        writer.writerow(_comparison_row())
    eval_file.write_text(json.dumps(_eval_row(), sort_keys=True) + "\n", encoding="utf-8")

    result = run_evaluation(comparison, eval_file, output_csv, report)

    assert result["rows"] == 4
    with output_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert REQUIRED_EXTENDED_COLUMNS <= set(rows[0])
    assert all(math.isfinite(float(row["multi_objective_score"])) for row in rows)
    assert "4-way Comparison Summary" in report.read_text(encoding="utf-8")
