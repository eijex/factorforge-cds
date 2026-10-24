"""Tests for Equivalence-Preserving Constrained Local Codon Rewriter (Job 440)."""

from __future__ import annotations

import pytest

from factorforge.analysis.metrics import (
    STANDARD_GENETIC_CODE,
    detect_homopolymers,
    translate_dna,
)
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.engines.dp_v2_3 import DPV230Optimizer
from factorforge.engines.local_rewriter import (
    ConstrainedLocalRewriter,
    build_4dp_pause_classes,
    compare_rewriter_vs_exact_dp,
    run_exact_dp_with_budget_trace,
    verify_small_instance_ground_truth,
)
from factorforge.utils.exceptions import UnsatisfiableDesignError


def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _synthetic_codon_weights() -> dict[str, float]:
    """Deterministic codon weights spanning high/medium/low frequencies per amino acid."""
    weights: dict[str, float] = {}
    by_aa: dict[str, list[str]] = {}
    for codon, aa in sorted(STANDARD_GENETIC_CODE.items()):
        if aa == "*":
            continue
        by_aa.setdefault(aa, []).append(codon)
    for aa, syns in sorted(by_aa.items()):
        syns_sorted = sorted(syns)
        n = len(syns_sorted)
        for idx, c in enumerate(syns_sorted):
            if n == 1:
                weights[c] = 1.0
            else:
                # Spread weights between 0.15 and 1.00
                weights[c] = round(0.15 + 0.85 * (idx / float(n - 1)), 4)
    return weights


def test_build_4dp_pause_classes_shifted_windows_and_validation() -> None:
    """Verifies W_k(L, delta) = [b_k + L - delta, b_k + L + delta] and input guards."""
    classes = build_4dp_pause_classes(
        boundaries=[12, 30],
        tunnel_length=35,
        half_width=2,
        weight=2.0,
        target_tau=0.45,
        protein_length=70,
    )
    assert "CUSTOM" in classes
    # b_1 = 12 -> 12 + 35 = 47 -> [45, 46, 47, 48, 49]
    # b_2 = 30 -> 30 + 35 = 65 -> [63, 64, 65, 66, 67]
    expected_positions = [45, 46, 47, 48, 49, 63, 64, 65, 66, 67]
    assert classes["CUSTOM"]["positions"] == expected_positions
    assert classes["CUSTOM"]["tunnel_length_assumption"] == 35
    assert classes["CUSTOM"]["half_width"] == 2

    # Clipping near protein end (1-based inclusive: 20 + 35 = 55 -> [53, 54, 55])
    clipped = build_4dp_pause_classes(
        boundaries=[20],
        tunnel_length=35,
        half_width=2,
        protein_length=55,
    )
    assert clipped["CUSTOM"]["positions"] == [53, 54, 55]

    with pytest.raises(ValueError, match="tunnel_length"):
        build_4dp_pause_classes([10], tunnel_length=-1)
    with pytest.raises(ValueError, match="half_width"):
        build_4dp_pause_classes([10], half_width=-1)
    with pytest.raises(ValueError, match="target_tau"):
        build_4dp_pause_classes([10], target_tau=1.5)


def test_rewriter_preserves_exact_translation_and_full_construct_constraints() -> None:
    """Every accepted rewrite preserves T(c') == T(c) == a and all sequence-level constraints."""
    weights = _synthetic_codon_weights()
    protein = "M" + "KELVAGSTNQRIFYD" * 3  # 46 aa
    dp = DPV230Optimizer(ramp_length_codons=15)
    rewriter = ConstrainedLocalRewriter(dp_optimizer=dp)

    res = rewriter.rewrite(
        protein_sequence=protein,
        codon_weights=weights,
        seed_mode="low_fitness_feasible",
        window_sizes=(1, 2, 3),
        include_paired_gc_swaps=True,
        target_gc_min=0.42,
        target_gc_max=0.52,
        upstream_context="GGGAAA",
        stop_codon="TAA",
        downstream_context="CCCGGG",
    )

    assert res["is_feasible"] is True
    assert res["translation_equivalent"] is True
    assert translate_dna(res["sequence"]) == protein
    assert res["biological_validation"] == "NOT_AVAILABLE"
    assert res["score"] >= res["initial_feasible_score"] - 1e-12
    assert res["rewrite_metrics"]["accepted_rewrites"] >= 1

    # Verify full construct has zero forbidden motifs or homopolymers
    full_construct = "GGGAAA" + res["sequence"] + "TAA" + "CCCGGG"
    for motif in get_canonical_forbidden_motifs():
        assert motif not in full_construct
        assert _reverse_complement(motif) not in full_construct
    assert detect_homopolymers(full_construct, max_run=dp.homopolymer_max_run + 1) == []

    # Verify every step in trace was monotonic and synonymous
    prev_score = res["initial_feasible_score"]
    for step in res["rewrite_metrics"]["trace"]:
        assert step["score_before"] == pytest.approx(prev_score, abs=1e-9)
        assert step["score_after"] >= step["score_before"] - 1e-12
        for pos, old_c, new_c in zip(step["positions"], step["old_codons"], step["new_codons"]):
            assert STANDARD_GENETIC_CODE[old_c] == protein[pos]
            assert STANDARD_GENETIC_CODE[new_c] == protein[pos]
        prev_score = step["score_after"]


def test_rewriter_rejects_non_synonymous_initial_seed() -> None:
    """A seed whose translated amino acids differ from protein_sequence must raise ValueError."""
    weights = _synthetic_codon_weights()
    rewriter = ConstrainedLocalRewriter()
    protein = "MKELVA"
    bad_dna = "ATG" + "TTT" * 5  # Translates to MFFFFF != MKELVA

    with pytest.raises(ValueError, match="amino acid translation equivalence"):
        rewriter.rewrite(
            protein_sequence=protein,
            codon_weights=weights,
            initial_coding_dna=bad_dna,
        )


def test_rewriter_repairs_infeasible_synonymous_seed_or_raises_when_required() -> None:
    """Verifies local seed repair when an initial synonymous sequence violates motifs/GC."""
    weights = _synthetic_codon_weights()
    protein = "MKELVAGSTNQRIFYD"
    dp = DPV230Optimizer(ramp_length_codons=6)
    rewriter = ConstrainedLocalRewriter(dp_optimizer=dp)

    # Construct a synonymous seed with very low GC (violates [0.45, 0.55])
    low_gc_codons = []
    for aa in protein:
        syns = sorted(
            [c for c, a in STANDARD_GENETIC_CODE.items() if a == aa],
            key=lambda c: (c.count("G") + c.count("C"), c),
        )
        low_gc_codons.append(syns[0])
    infeasible_seed = "".join(low_gc_codons)

    audit = rewriter.evaluate_candidate(
        protein_sequence=protein,
        coding_dna=infeasible_seed,
        codon_weights=weights,
        target_gc_min=0.45,
        target_gc_max=0.55,
    )
    assert audit["is_feasible"] is False

    with pytest.raises(UnsatisfiableDesignError, match="require_feasible_seed=True"):
        rewriter.rewrite(
            protein_sequence=protein,
            codon_weights=weights,
            initial_coding_dna=infeasible_seed,
            require_feasible_seed=True,
            target_gc_min=0.45,
            target_gc_max=0.55,
        )

    repaired = rewriter.rewrite(
        protein_sequence=protein,
        codon_weights=weights,
        initial_coding_dna=infeasible_seed,
        require_feasible_seed=False,
        target_gc_min=0.45,
        target_gc_max=0.55,
    )
    assert repaired["is_feasible"] is True
    assert repaired["rewrite_metrics"]["seed_repaired"] is True
    assert translate_dna(repaired["sequence"]) == protein


def test_coupled_gc_swap_escapes_single_site_gc_boundary_trap() -> None:
    """When global GC is locked at an exact integer count, W=1 cannot swap GC, but coupled swap can."""
    protein = "MEEEE"  # M (ATG: 1 GC) + 4 Glutamates (GAA: 1 GC, GAG: 2 GC)
    # Lock exact GC count = 7 out of 15 nt (7/15 = 0.4667):
    # ATG (1 GC) + two GAG (2 GC) + two GAA (1 GC) = 7 GC.
    # Because GAA (1 GC) and GAG (2 GC) have different GC counts, NO W=1 move can preserve GC=7.
    weights = _synthetic_codon_weights()
    weights["GAA"] = 0.20
    weights["GAG"] = 1.00

    # 1-based position 2 (0-indexed pos 1) has a pause favoring target_tau=0.20 (GAA),
    # while 0-indexed pos 4 (1-based pos 5) is a non-pause body codon favoring high weight (GAG).
    pause_classes = {
        "BOUNDARY": {
            "positions": [2],
            "weight": 4.0,
            "target_tau": 0.20,
        }
    }

    # Seed has GAG at pos 1 (bad for pause) and GAA at pos 4 (bad for body), total GC = 7:
    # pos 0: ATG (1 GC), pos 1: GAG (2 GC), pos 2: GAA (1 GC), pos 3: GAG (2 GC), pos 4: GAA (1 GC)
    seed_dna = "ATGGAGGAAGAGGAA"

    dp_tight = DPV230Optimizer(pause_classes=pause_classes, ramp_length_codons=1)
    rewriter_tight = ConstrainedLocalRewriter(dp_optimizer=dp_tight)

    # Global GC band [7/15, 7/15] -> min_gc_count = 7, max_gc_count = 7
    res_w1_only = rewriter_tight.rewrite(
        protein_sequence=protein,
        codon_weights=weights,
        initial_coding_dna=seed_dna,
        window_sizes=(1,),
        include_paired_gc_swaps=False,
        target_gc_min=7.0 / 15.0 - 0.001,
        target_gc_max=7.0 / 15.0 + 0.001,
    )
    res_with_coupled = rewriter_tight.rewrite(
        protein_sequence=protein,
        codon_weights=weights,
        initial_coding_dna=seed_dna,
        window_sizes=(1,),
        include_paired_gc_swaps=True,
        target_gc_min=7.0 / 15.0 - 0.001,
        target_gc_max=7.0 / 15.0 + 0.001,
    )

    assert res_w1_only["rewrite_metrics"]["accepted_rewrites"] == 0
    assert res_with_coupled["score"] > res_w1_only["score"] + 1e-6
    move_types = [m["move_type"] for m in res_with_coupled["rewrite_metrics"]["trace"]]
    assert "coupled_gc_swap" in move_types


def test_small_instance_exhaustive_ground_truth_verification() -> None:
    """Verifies that Exact DP matches brute-force ground truth and Rewriter matches on N=6."""
    weights = _synthetic_codon_weights()
    protein = "MKELVA"  # 1 * 2 * 2 * 6 * 4 * 4 = 384 synonymous candidates
    dp = DPV230Optimizer(ramp_length_codons=3)

    gt = verify_small_instance_ground_truth(
        protein_sequence=protein,
        codon_weights=weights,
        dp_optimizer=dp,
        target_gc_min=0.40,
        target_gc_max=0.56,
    )

    assert gt["total_synonymous_candidates"] == 384
    assert gt["feasible_candidates"] > 0
    assert gt["dp_matches_exhaustive_score"] is True
    assert gt["dp_matches_exhaustive_sequence"] is True
    assert gt["rewriter_feasible"] is True
    assert gt["rewriter_matches_exhaustive_score"] is True


def test_compare_rewriter_vs_exact_dp_across_1dp_3dp_and_4dp_modes() -> None:
    """Runs matched-objective, matched-constraint comparison across 1DP, 3DP, and 4DP (L=35)."""
    weights = _synthetic_codon_weights()
    protein = "M" + "KELVAGSTNQRIFYDCWPH" * 3  # 58 aa

    # 1. 1DP (no pause classes)
    dp_1d = DPV230Optimizer(pause_classes=None, ramp_length_codons=15)
    comp_1d = compare_rewriter_vs_exact_dp(
        protein_sequence=protein,
        codon_weights=weights,
        dp_optimizer=dp_1d,
        budgets=(100, 400, None),
        target_gc_min=0.42,
        target_gc_max=0.52,
    )
    assert comp_1d["biological_validation"] == "NOT_AVAILABLE"
    assert comp_1d["exact_dp"]["is_feasible"] is True
    scores_1d = [r["score"] for r in comp_1d["rewriter_runs"]]
    assert scores_1d[0] <= scores_1d[1] + 1e-12 <= scores_1d[2] + 1e-12
    assert scores_1d[-1] <= comp_1d["exact_dp"]["score"] + 1e-9

    # 2. 3DP (static boundary positions)
    pause_3d = {
        "BOUNDARY": {
            "positions": [18, 19, 20, 36, 37, 38],
            "weight": 2.0,
            "target_tau": 0.45,
        }
    }
    dp_3d = DPV230Optimizer(pause_classes=pause_3d, ramp_length_codons=15)
    comp_3d = compare_rewriter_vs_exact_dp(
        protein_sequence=protein,
        codon_weights=weights,
        dp_optimizer=dp_3d,
        budgets=(150, 600, None),
        target_gc_min=0.42,
        target_gc_max=0.52,
    )
    assert comp_3d["exact_dp"]["is_feasible"] is True
    for run in comp_3d["rewriter_runs"]:
        assert run["is_feasible"] is True
        assert run["translation_equivalent"] is True
        assert run["monotonic_trace_verified"] is True
        if run["budget_limit"] is not None:
            assert run["step_evaluations"] <= run["budget_limit"]

    # 3. 4DP (exit-tunnel shifted boundary windows W_k(L=35, delta=2))
    pause_4d = build_4dp_pause_classes(
        boundaries=[12, 20],
        tunnel_length=35,
        half_width=2,
        weight=2.0,
        target_tau=0.45,
        protein_length=len(protein),
    )
    dp_4d = DPV230Optimizer(pause_classes=pause_4d, ramp_length_codons=15)
    comp_4d = compare_rewriter_vs_exact_dp(
        protein_sequence=protein,
        codon_weights=weights,
        dp_optimizer=dp_4d,
        budgets=(150, 600, None),
        target_gc_min=0.42,
        target_gc_max=0.52,
    )
    assert comp_4d["exact_dp"]["is_feasible"] is True
    final_4d = comp_4d["rewriter_runs"][-1]
    assert final_4d["is_feasible"] is True
    assert final_4d["score_improvement_over_seed"] > 0.0
    assert final_4d["optimality_gap_vs_dp"] >= -1e-9


def test_run_exact_dp_with_budget_trace_reports_positive_step_evaluations() -> None:
    """Verifies exact DP transition step accounting for matched-budget comparisons."""
    weights = _synthetic_codon_weights()
    protein = "MKELVAGSTNQRIFYD"
    dp = DPV230Optimizer(ramp_length_codons=5)
    summary = run_exact_dp_with_budget_trace(
        dp_optimizer=dp,
        protein_sequence=protein,
        codon_weights=weights,
        target_gc_min=0.40,
        target_gc_max=0.55,
    )
    assert summary["is_feasible"] is True
    assert summary["translation_equivalent"] is True
    assert summary["step_evaluations"] > len(protein)
    assert summary["peak_active_states"] >= 1
