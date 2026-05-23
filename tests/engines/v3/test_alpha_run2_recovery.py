"""Tests for Run 4B recovery data prep and smoke evaluation."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest
yaml = pytest.importorskip("yaml", reason="pyyaml not installed — skipping ML tests")

ROOT = Path(__file__).resolve().parents[3]
DATA_PREP = ROOT / "scripts" / "1_data_preparation"
EVALUATION = ROOT / "scripts" / "4_evaluation"
for path in [DATA_PREP, EVALUATION]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from generate_run4b_mixed_pseudolabels import (  # noqa: E402
    AVAILABLE_V2_PROFILES,
    build_mixed_teacher_row,
    generate_run4b_mixed_pseudolabels,
)
from run4b_smoke_compare import summarize_result  # noqa: E402


RUN4B_SMOKE_CONFIGS = [
    ROOT / "configs" / "v3_training_config_run4b_gc1.yml",
    ROOT / "configs" / "v3_training_config_run4b_gc2_low2.yml",
    ROOT / "configs" / "v3_training_config_run4b_gc2_low2_cai010.yml",
]
RUN4B_FULL_CONFIG_B = ROOT / "configs" / "v3_training_config_alpha_run2.yml"


def test_v2_profile_availability_contract() -> None:
    assert "gc_target" in AVAILABLE_V2_PROFILES
    assert "local_gc_balanced" not in AVAILABLE_V2_PROFILES
    assert "motif_clean" not in AVAILABLE_V2_PROFILES


def test_mixed_pseudolabel_generates_valid_cds_for_each_amino_acid(tmp_path: Path) -> None:
    fasta = tmp_path / "proteins.fasta"
    output = tmp_path / "run4b.jsonl"
    fasta.write_text(">all_aa\nMACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8")

    counts = generate_run4b_mixed_pseudolabels(fasta, output)
    row = json.loads(output.read_text(encoding="utf-8").strip())

    assert counts["written"] == 1
    assert row["amino_acid_sequence"] == "MACDEFGHIKLMNPQRSTVWY"
    assert len(row["dna_sequence"]) == len(row["amino_acid_sequence"]) * 3
    assert row["validator_result"]["passed"] is True
    assert row["selected_teacher_type"] in row["candidate_metrics"]
    assert isinstance(row["selection_score"], float)


def test_selected_teacher_records_type_and_metrics() -> None:
    row = build_mixed_teacher_row("low_gc_prone", "MKIYNNK")

    assert row["selected_teacher_type"]
    assert "high_cai" in row["candidate_metrics"]
    assert "gc_target" in row["candidate_metrics"]
    assert "gc" in row["candidate_metrics"][row["selected_teacher_type"]]
    assert row["validator_result"]["amino_acid_identity"] == 1.0


def test_gc_target_teacher_has_higher_gc_than_high_cai_when_feasible() -> None:
    row = build_mixed_teacher_row("low_gc_prone", "MKIYNNK")
    high_cai_gc = row["candidate_metrics"]["high_cai"]["gc"]
    gc_target_gc = row["candidate_metrics"]["gc_target"]["gc"]

    assert gc_target_gc > high_cai_gc


def test_run4b_smoke_configs_are_bounded_smoke_configs() -> None:
    for path in RUN4B_SMOKE_CONFIGS:
        cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert cfg["smoke"] is True
        assert cfg["training"]["max_steps"] <= 2000
        assert cfg["loss"]["gc_weight"] in {1.0, 2.0}
        assert cfg["loss"]["expected_log_cai_weight"] in {0.05, 0.10}


def test_run4b_full_config_b_contract() -> None:
    cfg = yaml.safe_load(RUN4B_FULL_CONFIG_B.read_text(encoding="utf-8"))

    assert cfg["smoke"] is False
    assert cfg["training"]["max_steps"] == 20000
    assert cfg["training"]["checkpoint_every"] == 5000
    assert cfg["loss"]["gc_weight"] == 2.0
    assert cfg["loss"]["gc_lambda_low"] == 1.0
    assert cfg["loss"]["expected_log_cai_weight"] == 0.05
    assert cfg["data"]["train_file"] == "data/training/run4b_pseudolabels_train.jsonl"
    assert cfg["data"]["eval_file"] == "data/training/run4b_pseudolabels_eval.jsonl"
    assert cfg["data"]["pseudo_label_path"] == cfg["data"]["train_file"]


def test_run4b_evaluation_parses_outputs(tmp_path: Path) -> None:
    comparison = tmp_path / "run4_comparison.csv"
    loss = tmp_path / "run4_loss_log.csv"
    summary = tmp_path / "run4_summary.json"

    with comparison.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "protein_id",
                "v2_cai",
                "v3_cai",
                "v2_gc_global",
                "v3_gc_global",
                "v3_amino_acid_identity",
                "v3_validator_pass",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "protein_id": "p1",
                "v2_cai": "0.8",
                "v3_cai": "0.79",
                "v2_gc_global": "35.0",
                "v3_gc_global": "42.0",
                "v3_amino_acid_identity": "1.0",
                "v3_validator_pass": "True",
            }
        )
        writer.writerow(
            {
                "protein_id": "p2",
                "v2_cai": "0.7",
                "v3_cai": "0.69",
                "v2_gc_global": "34.0",
                "v3_gc_global": "44.0",
                "v3_amino_acid_identity": "1.0",
                "v3_validator_pass": "True",
            }
        )
    with loss.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "step",
                "split",
                "expected_GC",
                "bounded_GC_penalty",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "step": "2000",
                "split": "eval",
                "expected_GC": "0.42",
                "bounded_GC_penalty": "0.0",
            }
        )
    summary.write_text("{}", encoding="utf-8")

    result = summarize_result(tmp_path)

    assert result["decoded_count"] == 2
    assert result["validator_pass_rate"] == 1.0
    assert result["v3_global_gc_mean"] == 43.0
    assert result["expected_gc_mean"] == 0.42
    assert result["verdict"] == "pass"
