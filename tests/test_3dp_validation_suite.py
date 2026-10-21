"""Tests for FactorForge 3DP Validation Suite.

Verifies:
1. Lambda=0 Equivalence: 3DP collapses to exact 1DP.
2. Hard Constraint Regression: Translation identity, GC band, Type IIS motifs preserved.
3. Cost-of-Conditioning & Structural Gain calculation.
4. Edit Localization metric.
5. Shuffled Annotation negative control (causal intervention shift).
6. Boundary Perturbation stability.
7. Pareto Lambda Sweep profiling.
"""

import pytest

from factorforge.validation.tradeoff_suite import (
    check_lambda_zero_equivalence,
    compute_tradeoff_comparison,
    run_pareto_lambda_sweep,
    run_shuffled_structure_control,
    run_boundary_perturbation_test,
)


@pytest.fixture
def benthamiana_weights():
    # FactorForge engineering tests must run without the private AgentOS package.
    from factorforge.analysis.metrics import load_codon_usage_table

    return load_codon_usage_table().codon_weights


@pytest.fixture
def test_antibody_fragment():
    # 30 amino acid synthetic antibody segment
    # Positions 11-18 represent a structural linker
    return (
        "DIQMTQSPSS"  # 1-10
        "ASTLQSGV"    # 11-18 (Linker)
        "PSRFSGSGSG"  # 19-28
        "TD"          # 29-30
    )


def test_suite_lambda_zero_equivalence(test_antibody_fragment, benthamiana_weights):
    """Test 1: λ=0 Equivalence Check."""
    pause_pos = list(range(11, 19))
    is_equivalent = check_lambda_zero_equivalence(
        protein=test_antibody_fragment,
        pause_positions=pause_pos,
        codon_weights=benthamiana_weights,
        target_gc_min=0.40,
        target_gc_max=0.48,
        forbidden_motifs=["GGTCTC"],
        target_adaptiveness=0.40,
    )
    assert is_equivalent is True


def test_suite_hard_constraints_and_tradeoff_quantification(test_antibody_fragment, benthamiana_weights):
    """Tests 2, 3, 4, 5: Hard constraints, trade-off metrics, and localization."""
    pause_pos = list(range(11, 19))
    comp = compute_tradeoff_comparison(
        protein=test_antibody_fragment,
        pause_positions=pause_pos,
        codon_weights=benthamiana_weights,
        target_gc_min=0.40,
        target_gc_max=0.48,
        forbidden_motifs=["GGTCTC", "GAGACC"],
        lambda_weight=2.0,
        target_adaptiveness=0.40,
    )

    # Test 2: Hard constraints must pass
    assert comp.hard_constraints_pass is True

    # Test 3 & 4: Cost of conditioning and structural gain
    # 3DP should sacrifice some 1D CAI (Δ1D <= 0) to gain structural fitness (Δ3D > 0)
    assert comp.delta_3d_gain > 0.0, "3DP must achieve positive structural context gain"
    assert comp.delta_1d_cai <= 0.001, "Cost of conditioning: CAI should reflect trade-off"

    # Test 5: Localization
    assert comp.total_edits > 0
    assert comp.target_zone_edits > 0
    assert comp.localization_ratio >= 0.50, f"Edits should localize to target zone: {comp.localization_ratio}"

    # Verify dashboard table format
    assert "Structural Conditioning Cost:" in comp.dashboard_table
    assert "Structural Context Gain:" in comp.dashboard_table
    assert "Localization Ratio:" in comp.dashboard_table
    print("\n" + comp.dashboard_table)


def test_suite_shuffled_annotation_negative_control(test_antibody_fragment, benthamiana_weights):
    """Test 6: Shuffled-Structure Negative Control."""
    nominal_pos = list(range(11, 19))  # Positions 11-18
    shuffled_pos = list(range(2, 10))   # Positions 2-9 (Shuffled)

    res = run_shuffled_structure_control(
        protein=test_antibody_fragment,
        nominal_positions=nominal_pos,
        shuffled_positions=shuffled_pos,
        codon_weights=benthamiana_weights,
        target_gc_min=0.40,
        target_gc_max=0.48,
        lambda_weight=2.0,
        target_adaptiveness=0.40,
    )

    assert res.causal_shift_confirmed is True
    assert res.shuffled_edits_in_shuffled_zone > 0


def test_suite_boundary_perturbation_stability(test_antibody_fragment, benthamiana_weights):
    """Test 7: Boundary Perturbation Sensitivity."""
    nominal_pos = list(range(11, 19))
    res = run_boundary_perturbation_test(
        protein=test_antibody_fragment,
        nominal_positions=nominal_pos,
        offsets=(-1, 0, 1),
        codon_weights=benthamiana_weights,
        target_gc_min=0.40,
        target_gc_max=0.48,
        lambda_weight=2.0,
        target_adaptiveness=0.40,
    )
    assert res.is_stable is True
    for offset, data in res.perturbations.items():
        assert data["structure_fit"] > 0.0


def test_suite_pareto_lambda_sweep(test_antibody_fragment, benthamiana_weights):
    """Test 8: Pareto λ Sweep."""
    pause_pos = list(range(11, 19))
    sweep = run_pareto_lambda_sweep(
        protein=test_antibody_fragment,
        pause_positions=pause_pos,
        codon_weights=benthamiana_weights,
        lambda_values=[0.0, 0.5, 1.0, 2.0],
        target_gc_min=0.40,
        target_gc_max=0.48,
        target_adaptiveness=0.40,
    )
    assert len(sweep) == 4
    # At lambda = 0.0: edits should be 0
    assert sweep[0].lambda_val == 0.0
    assert sweep[0].total_edits_from_1dp == 0

    # As lambda increases, structural fit should increase or stay high
    assert sweep[-1].structure_fit >= sweep[0].structure_fit
