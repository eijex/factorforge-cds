"""Unit and Regression Tests for Matched-Condition DP vs 3DP vs 4DP Comparison (Job 508 / Job 509)."""

from __future__ import annotations

import pytest

from factorforge.analysis.metrics import load_codon_usage_table, translate_dna
from factorforge.validation.dp_4dp_comparison import (
    BENCHMARK_PROTEIN_PANEL,
    make_not_run_invariants_record,
    run_dp_4dp_benchmark_panel,
    run_matched_dp_comparison,
    verify_dp_4dp_collapse_invariants,
    verify_small_instance_exhaustive_dp_4dp,
)


@pytest.fixture(scope="module")
def default_weights() -> dict[str, float]:
    return load_codon_usage_table().codon_weights


def test_small_instance_exhaustive_ground_truth_for_1dp_3dp_4dp(
    default_weights: dict[str, float],
) -> None:
    """Verifies exact global optimality of 1DP, 3DP, and 4DP against brute-force enumeration on N=8."""
    gt = verify_small_instance_exhaustive_dp_4dp(
        protein_sequence="MVSKGEEL",
        boundaries=[2],
        codon_weights=default_weights,
        tunnel_length=4,
        half_width=1,
        weight=2.5,
        target_tau=0.45,
    )

    assert gt["execution_status"] == "EXECUTED"
    assert gt["static_3dp_window_positions"] == [1, 2, 3]
    assert gt["shifted_4dp_window_positions"] == [5, 6, 7]
    assert gt["all_modes_match_exhaustive_global_optimum_score"] is True
    assert gt["all_modes_match_exhaustive_global_optimum_sequence"] is True

    for mode_key in ("1dp", "3dp", "4dp"):
        m = gt["modes"][mode_key]
        assert m["feasible_candidates"] > 0
        assert m["dp_matches_exhaustive_score"] is True
        assert m["dp_matches_exhaustive_sequence"] is True


def test_matched_dp_comparison_hard_constraints_and_candidate_self_preference(
    default_weights: dict[str, float],
) -> None:
    """Verifies hard constraints, 3-candidate self-preference, independent state replay, and memory scope."""
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
        independent_state_replay=True,
    )

    assert res["evidence_classification"] == "UNVERIFIED_PUBLIC_CORPUS"
    assert res["biological_validation"] == "NOT_AVAILABLE"
    assert res["all_constraints_satisfied"] is True
    assert res["candidate_self_preference_verified"] is True
    assert (
        res["cross_objective_matrix_scope"]
        == "THREE_CANDIDATE_SELF_PREFERENCE_ONLY_NOT_EXHAUSTIVE_GLOBAL_OPTIMALITY_PROOF"
    )
    assert res["state_graph_parity_verified"] is True
    assert res["state_metrics_measurement_mode"] == "INDEPENDENT_STATE_GRAPH_REPLAY"

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
        assert m["python_tracemalloc_peak_kb"] > 0.0
        assert (
            m["memory_measurement_scope"]
            == "PYTHON_TRACEMALLOC_HEAP_ONLY_EXCLUDES_NATIVE_RSS"
        )
        assert m["step_evaluations"] > 0
        assert m["peak_active_states"] >= 1
        assert m["state_metrics_measurement_mode"] == "INDEPENDENT_STATE_GRAPH_REPLAY"

    # Verify independent state replay produces identical counts across all 3 modes
    assert (
        res["modes"]["1dp"]["step_evaluations"]
        == res["modes"]["3dp"]["step_evaluations"]
        == res["modes"]["4dp"]["step_evaluations"]
    )

    # Verify 3x3 cross-objective matrix self-preference among the 3 generated candidates
    matrix = res["cross_objective_matrix"]
    for obj_mode in ("1dp", "3dp", "4dp"):
        own_val = matrix[obj_mode][obj_mode]
        for cand_mode in ("1dp", "3dp", "4dp"):
            assert own_val >= matrix[obj_mode][cand_mode] - 1e-6


def test_shared_state_replay_mode_explicitly_labeled(
    default_weights: dict[str, float],
) -> None:
    """Verifies that when independent_state_replay=False, modes are explicitly labeled as shared reuse."""
    fixture = BENCHMARK_PROTEIN_PANEL[4]
    res = run_matched_dp_comparison(
        protein_sequence=fixture["protein"],
        boundaries=fixture["boundaries"],
        codon_weights=default_weights,
        independent_state_replay=False,
    )
    assert res["state_metrics_measurement_mode"] == "SHARED_COMMON_GRAPH_REUSE"
    assert res["modes"]["1dp"]["state_metrics_measurement_mode"] == "INDEPENDENT_STATE_GRAPH_REPLAY"
    assert res["modes"]["3dp"]["state_metrics_measurement_mode"] == "SHARED_COMMON_GRAPH_REUSE"
    assert res["modes"]["4dp"]["state_metrics_measurement_mode"] == "SHARED_COMMON_GRAPH_REUSE"


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

    assert inv["execution_status"] == "EXECUTED"
    assert inv["lambda_zero_collapse_verified"] is True
    assert inv["tunnel_zero_collapse_verified"] is True
    assert inv["out_of_bounds_exclusion_verified"] is True
    assert inv["all_invariants_verified"] is True
    assert inv["details"]["oob_shifted_4dp_positions"] == []


def test_not_run_invariants_never_claim_true() -> None:
    """Verifies that skipped invariant checks return NOT_RUN and None rather than True."""
    rec = make_not_run_invariants_record("Skipped for runtime")
    assert rec["execution_status"] == "NOT_RUN"
    assert rec["lambda_zero_collapse_verified"] is None
    assert rec["tunnel_zero_collapse_verified"] is None
    assert rec["out_of_bounds_exclusion_verified"] is None
    assert rec["all_invariants_verified"] is None


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
        independent_state_replay=False,
    )

    setup = res["problem_setup"]
    assert setup["static_3dp_window_positions"] == [26, 27, 28, 29, 30]
    assert setup["shifted_4dp_window_positions"] == [61, 62, 63, 64, 65]

    m_1dp = res["modes"]["1dp"]
    m_3dp = res["modes"]["3dp"]
    m_4dp = res["modes"]["4dp"]

    assert m_3dp["mean_abs_tau_dev_in_3dp_window"] <= m_1dp["mean_abs_tau_dev_in_3dp_window"] + 1e-6
    assert m_3dp["mean_abs_tau_dev_in_3dp_window"] <= m_4dp["mean_abs_tau_dev_in_3dp_window"] + 1e-6

    assert m_4dp["mean_abs_tau_dev_in_4dp_window"] <= m_1dp["mean_abs_tau_dev_in_4dp_window"] + 1e-6
    assert m_4dp["mean_abs_tau_dev_in_4dp_window"] <= m_3dp["mean_abs_tau_dev_in_4dp_window"] + 1e-6

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
        independent_state_replay=False,
    )

    assert res["problem_setup"]["static_3dp_window_positions"] == [28, 29, 30, 31, 32]
    assert res["problem_setup"]["shifted_4dp_window_positions"] == []
    assert res["modes"]["4dp"]["sequence"] == res["modes"]["1dp"]["sequence"]
    assert res["pairwise_candidate_divergence"]["1dp_vs_4dp"]["codon_hamming_distance"] == 0
    assert res["pairwise_candidate_divergence"]["1dp_vs_3dp"]["codon_hamming_distance"] > 0


def test_full_benchmark_panel_honest_invariant_accounting(
    default_weights: dict[str, float],
) -> None:
    """Runs the 5-protein benchmark panel with verify_invariants_for_all=False and checks NOT_RUN accounting."""
    panel_report = run_dp_4dp_benchmark_panel(
        codon_weights=default_weights,
        verify_invariants_for_all=False,
        independent_state_replay=True,
    )
    assert panel_report["benchmark_target_count"] == 5
    assert panel_report["all_targets_feasible"] is True
    assert panel_report["all_targets_candidate_self_preference_verified"] is True
    assert (
        panel_report["small_instance_exhaustive_ground_truth"][
            "all_modes_match_exhaustive_global_optimum_score"
        ]
        is True
    )
    # 3 short fixtures (<= 100 aa) executed, 2 longer fixtures (> 100 aa) marked NOT_RUN (None)
    assert panel_report["invariants_executed_target_count"] == 3
    assert panel_report["invariants_not_run_target_count"] == 2
    assert panel_report["all_executed_invariants_verified"] is True
    assert panel_report["all_targets_invariants_verified"] is None

    for entry in panel_report["entries"]:
        assert "fixture_kind" in entry
        assert "sequence_provenance" in entry
        assert "boundary_provenance" in entry
        length_aa = entry["comparison"]["problem_setup"]["protein_length_aa"]
        if length_aa <= 100:
            assert entry["invariants"]["execution_status"] == "EXECUTED"
            assert entry["invariants"]["all_invariants_verified"] is True
        else:
            assert entry["invariants"]["execution_status"] == "NOT_RUN"
            assert entry["invariants"]["all_invariants_verified"] is None
