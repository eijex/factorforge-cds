"""Unit and Regression Tests for Matched-Condition DP vs 3DP vs 4DP Comparison (Job 508)."""

from __future__ import annotations

import pytest

from factorforge.analysis.metrics import load_codon_usage_table, translate_dna
from factorforge.validation.dp_4dp_comparison import (
    BENCHMARK_PROTEIN_PANEL,
    run_dp_4dp_benchmark_panel,
    run_matched_dp_comparison,
    verify_dp_4dp_collapse_invariants,
)


@pytest.fixture(scope="module")
def default_weights() -> dict[str, float]:
    return load_codon_usage_table().codon_weights


def test_matched_dp_comparison_hard_constraints_and_diagonal_optimality(
    default_weights: dict[str, float],
) -> None:
    """Verifies that 1DP, 3DP, and 4DP satisfy all hard constraints and achieve diagonal global optimality."""
    protein = BENCHMARK_PROTEIN_PANEL[0]["protein"]
    boundaries = BENCHMARK_PROTEIN_PANEL[0]["boundaries"]

    res = run_matched_dp_comparison(
        protein_sequence=protein,
        boundaries=boundaries,
        codon_weights=default_weights,
        tunnel_length=35,
        half_width=2,
        weight=2.5,
        target_tau=0.45,
        upstream_context="CACC",
        stop_codon="TAA",
        downstream_context="GCTT",
    )

    assert res["evidence_classification"] == "UNVERIFIED_PUBLIC_CORPUS"
    assert res["biological_validation"] == "NOT_AVAILABLE"
    assert res["all_constraints_satisfied"] is True
    assert res["all_diagonal_optimal"] is True

    for mode_key in ("1dp", "3dp", "4dp"):
        m = res["modes"][mode_key]
        assert m["is_feasible"] is True
        assert m["translation_equivalent"] is True
        assert translate_dna(m["sequence"]) == protein
        assert m["automaton_valid"] is True
        assert m["initiation_gc_feasible"] is True
        assert m["global_gc_feasible"] is True
        assert m["homopolymers_detected"] == 0
        assert m["violations"] == []
        assert m["wall_time_ms"] > 0.0
        assert m["peak_memory_kb"] > 0.0
        assert m["step_evaluations"] > 0
        assert m["peak_active_states"] >= 1

    # Verify 3x3 cross-objective matrix diagonal dominance
    matrix = res["cross_objective_matrix"]
    for obj_mode in ("1dp", "3dp", "4dp"):
        own_val = matrix[obj_mode][obj_mode]
        for cand_mode in ("1dp", "3dp", "4dp"):
            assert own_val >= matrix[obj_mode][cand_mode] - 1e-6


def test_dp_4dp_mathematical_collapse_and_oob_invariants(
    default_weights: dict[str, float],
) -> None:
    """Verifies lambda=0 collapse, L=0 collapse (4DP == 3DP), and C-terminal out-of-bounds exclusion."""
    protein = BENCHMARK_PROTEIN_PANEL[0]["protein"]
    boundaries = BENCHMARK_PROTEIN_PANEL[0]["boundaries"]

    inv = verify_dp_4dp_collapse_invariants(
        protein_sequence=protein,
        boundaries=boundaries,
        codon_weights=default_weights,
        tunnel_length=35,
        half_width=2,
        target_tau=0.45,
    )

    assert inv["lambda_zero_collapse_verified"] is True
    assert inv["tunnel_zero_collapse_verified"] is True
    assert inv["out_of_bounds_exclusion_verified"] is True
    assert inv["all_invariants_verified"] is True
    assert inv["details"]["oob_shifted_4dp_positions"] == []


def test_spatial_shift_between_3dp_and_4dp_windows(
    default_weights: dict[str, float],
) -> None:
    """Verifies that 3DP and 4DP localize their synonymous attenuation to D_k vs W_k(L, delta) respectively."""
    protein = BENCHMARK_PROTEIN_PANEL[1]["protein"]
    boundaries = [28]

    res = run_matched_dp_comparison(
        protein_sequence=protein,
        boundaries=boundaries,
        codon_weights=default_weights,
        tunnel_length=35,
        half_width=2,
        weight=3.0,
        target_tau=0.40,
    )

    setup = res["problem_setup"]
    assert setup["static_3dp_window_positions"] == [26, 27, 28, 29, 30]
    assert setup["shifted_4dp_window_positions"] == [61, 62, 63, 64, 65]

    m_1dp = res["modes"]["1dp"]
    m_3dp = res["modes"]["3dp"]
    m_4dp = res["modes"]["4dp"]

    # In the static 3DP window [26..30], 3DP should achieve closer target_tau alignment than 1DP and 4DP
    assert m_3dp["mean_abs_tau_dev_in_3dp_window"] <= m_1dp["mean_abs_tau_dev_in_3dp_window"] + 1e-6
    assert m_3dp["mean_abs_tau_dev_in_3dp_window"] <= m_4dp["mean_abs_tau_dev_in_3dp_window"] + 1e-6

    # In the shifted 4DP window [61..65], 4DP should achieve closer target_tau alignment than 1DP and 3DP
    assert m_4dp["mean_abs_tau_dev_in_4dp_window"] <= m_1dp["mean_abs_tau_dev_in_4dp_window"] + 1e-6
    assert m_4dp["mean_abs_tau_dev_in_4dp_window"] <= m_3dp["mean_abs_tau_dev_in_4dp_window"] + 1e-6

    # Pairwise divergence between 3DP and 4DP should be non-zero when windows are disjoint
    diff_3_4 = res["pairwise_candidate_divergence"]["3dp_vs_4dp"]
    assert diff_3_4["codon_hamming_distance"] > 0
    assert (diff_3_4["edits_in_3dp_window"] + diff_3_4["edits_in_4dp_window"]) > 0


def test_cterminal_clipping_fixture_4dp_collapses_to_1dp_while_3dp_remains_active(
    default_weights: dict[str, float],
) -> None:
    """When b_k=30 on a 55 aa protein with L=35, delta=2, W_k=[63..67] > 55 is empty so 4DP == 1DP while 3DP != 1DP."""
    fixture = BENCHMARK_PROTEIN_PANEL[4]
    res = run_matched_dp_comparison(
        protein_sequence=fixture["protein"],
        boundaries=fixture["boundaries"],
        codon_weights=default_weights,
        tunnel_length=35,
        half_width=2,
        weight=2.5,
        target_tau=0.45,
    )

    assert res["problem_setup"]["static_3dp_window_positions"] == [28, 29, 30, 31, 32]
    assert res["problem_setup"]["shifted_4dp_window_positions"] == []
    assert res["modes"]["4dp"]["sequence"] == res["modes"]["1dp"]["sequence"]
    assert res["pairwise_candidate_divergence"]["1dp_vs_4dp"]["codon_hamming_distance"] == 0
    assert res["pairwise_candidate_divergence"]["1dp_vs_3dp"]["codon_hamming_distance"] > 0


def test_full_benchmark_panel_execution(default_weights: dict[str, float]) -> None:
    """Runs the 5-protein benchmark panel and verifies 100% feasibility, optimality, and invariant pass rates."""
    panel_report = run_dp_4dp_benchmark_panel(codon_weights=default_weights)
    assert panel_report["benchmark_target_count"] == 5
    assert panel_report["all_targets_feasible"] is True
    assert panel_report["all_targets_diagonal_optimal"] is True
    assert panel_report["all_targets_invariants_verified"] is True
