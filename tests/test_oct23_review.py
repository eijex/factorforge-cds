"""Synthetic regressions for Job441; no partner sequences or cloud artifacts."""

import json
from pathlib import Path

import pytest

from factorforge.engines import local_rewriter
from factorforge.rules.registry import RuleRegistry
from factorforge.sops.config import SopProfile


@pytest.mark.parametrize(
    "kwargs",
    [
        {"weight": True},
        {"target_tau": False},
        {"protein_length": -1},
        {"protein_length": True},
        {"class_name": ""},
    ],
)
def test_pause_window_invalid_types(kwargs):
    with pytest.raises(ValueError):
        local_rewriter.build_4dp_pause_classes([1], **kwargs)


def test_reference_constraint_mismatch_is_not_compared(monkeypatch):
    monkeypatch.setattr(
        local_rewriter, "run_exact_dp_with_budget_trace", lambda **kwargs: {"is_feasible": False}
    )
    with pytest.raises(ValueError, match="comparison rejected"):
        local_rewriter.compare_rewriter_vs_exact_dp("M", {"ATG": 1.0})


def test_actual_dp_ramp_relaxation_is_rejected_as_strict_comparison():
    optimizer = local_rewriter.DPV230Optimizer(
        ramp_length_codons=2,
        forbidden_motifs={"AAG", "CGT", "CGA", "CGG", "AGA", "AGG"},
        include_reverse_complement=False,
    )
    weights = {c: 1.0 for codons in local_rewriter.AA_TO_CODONS.values() for c in codons}
    with pytest.raises(ValueError, match="strict local constraint"):
        local_rewriter.compare_rewriter_vs_exact_dp(
            "KR", weights, dp_optimizer=optimizer, target_gc_min=0, target_gc_max=1
        )


def test_budget_comparison_declares_seed_and_estimate():
    result = local_rewriter.compare_rewriter_vs_exact_dp(
        "MA",
        {"ATG": 1.0, "GCT": 0.5, "GCC": 1.0, "GCA": 0.7, "GCG": 0.9},
        budgets=(30,),
        target_gc_min=0,
        target_gc_max=1,
    )
    assert result["matched_budget_comparison"] is False
    assert result["exact_dp"]["step_evaluations_kind"] == "SEPARATE_STATE_GRAPH_REPLAY_ESTIMATE"
    row = result["rewriter_runs"][0]
    assert (
        row["total_step_evaluations_including_seed"]
        == row["step_evaluations"] + row["seed_step_evaluations"]
    )
    assert row["step_evaluations"] <= 30


def test_collaboration_profile_is_draft_and_registry_loadable():
    path = Path(__file__).parents[1] / "configs/sops/plantform_benthamiana_production.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    profile = SopProfile(data, known_rule_ids=[r.rule_id for r in RuleRegistry().list_rules()])
    assert profile.status == "INTERNAL_REVIEW_DRAFT"
    assert data["approval_status"] == "NOT_DOCUMENTED_IN_LOCAL_REVIEW"
    assert data["authority_scope"] == "CONFIGURED_SEQUENCE_CHECKS_ONLY"
