import itertools
import math
import pytest
from factorforge.engines.dp_v2 import AA_TO_CODONS, DPV2Optimizer
from factorforge.engines.structural_pause import StructuralPauseDPOptimizer


def optimizer(**kwargs):
    return StructuralPauseDPOptimizer(
        positions={1},
        target_adaptiveness=0.4,
        penalty_weight=5.0,
        structure_reference="synthetic:manual:0.1.0",
        forbidden_motifs=[],
        **kwargs,
    )


def test_exact_objective_matches_exhaustive_search_and_cai_is_separate():
    protein = "KK"
    weights = {"AAA": 1.0, "AAG": 0.4}
    result = optimizer().optimize(protein, weights, 0, 1)
    scores = {
        "".join(c): sum(
            math.log(weights[x]) - (5 * abs(weights[x] - 0.4) if i == 0 else 0)
            for i, x in enumerate(c)
        )
        for c in itertools.product(AA_TO_CODONS["K"], repeat=2)
    }
    assert result["score"] == pytest.approx(max(scores.values()))
    assert result["sequence"] == "AAGAAA"
    assert result["cai"] == pytest.approx(math.sqrt(0.4))
    assert result["experimental"] is True


def test_empty_positions_reproduce_baseline():
    weights = {"AAA": 1.0, "AAG": 0.4}
    baseline = DPV2Optimizer(forbidden_motifs=[]).optimize("KK", weights, 0, 1)
    opt = StructuralPauseDPOptimizer(
        positions=[],
        target_adaptiveness=0.4,
        penalty_weight=10,
        structure_reference="synthetic",
        forbidden_motifs=[],
    )
    result = opt.optimize("KK", weights, 0, 1)
    for key in baseline:
        assert result[key] == baseline[key]


def test_gc_and_motif_constraints_remain_enforced():
    opt = StructuralPauseDPOptimizer(
        positions={1},
        target_adaptiveness=0.4,
        penalty_weight=5,
        structure_reference="synthetic",
        forbidden_motifs=["AAGAAA"],
    )
    result = opt.optimize("KK", {"AAA": 1, "AAG": 0.4}, 1 / 6, 1 / 6)
    assert result["sequence"] != "AAGAAA"
    assert result["gc_feasible"]


@pytest.mark.parametrize("positions", [{0}, {True}, {3}])
def test_invalid_positions_are_rejected(positions):
    with pytest.raises(ValueError):
        StructuralPauseDPOptimizer(
            positions=positions,
            target_adaptiveness=0.4,
            penalty_weight=5,
            structure_reference="synthetic",
        ).optimize("KK", {"AAA": 1})


def test_nonfinite_weights_rejected():
    with pytest.raises(ValueError, match="finite"):
        optimizer().optimize("KK", {"AAA": float("nan")})


def test_experimental_path_does_not_return_infeasible_fallback():
    with pytest.raises(RuntimeError, match="GC band"):
        optimizer().optimize("KK", {"AAA": 1, "AAG": 0.4}, 0.9, 1)
