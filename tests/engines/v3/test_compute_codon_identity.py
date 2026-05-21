"""Tests for v3/v2 gc_target codon identity calculation."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
EVALUATION = ROOT / "scripts" / "4_evaluation"
if str(EVALUATION) not in sys.path:
    sys.path.insert(0, str(EVALUATION))

from compute_codon_identity import (  # noqa: E402
    codon_identity_stats,
    run_identity,
)


def test_codon_identity_stats_known_pair() -> None:
    identity_pct, different_positions = codon_identity_stats("ATGGCCAAA", "ATGGCTAAA")

    assert identity_pct == pytest.approx(66.6666666667)
    assert different_positions == 1


def test_run_identity_writes_report_and_summary(tmp_path: Path) -> None:
    v3_decoded = tmp_path / "v3_decoded_eval.csv"
    output_csv = tmp_path / "codon_identity_report.csv"
    summary_json = tmp_path / "codon_identity_summary.json"

    with v3_decoded.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "protein_id",
                "amino_acid_sequence",
                "v3_decoded_cds",
                "validator_passed",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "protein_id": "p1",
                "amino_acid_sequence": "MA",
                "v3_decoded_cds": "ATGGCC",
                "validator_passed": "True",
            }
        )

    summary = run_identity(v3_decoded, output_csv, summary_json)

    assert summary["count"] == 1
    assert output_csv.exists()
    assert json.loads(summary_json.read_text(encoding="utf-8"))["count"] == 1
    with output_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["protein_id"] == "p1"
    assert rows[0]["v2_gc_target_cds"]
