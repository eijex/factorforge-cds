"""Regression checks for unsupported objective options and diagnostic claim scope."""

import pytest
from factorforge.engines.dp_v2_4 import DPV240Optimizer, TransferMatrixFeasibilityGuard
from factorforge.engines.sllm.automaton import AutomatonCompiler


@pytest.mark.parametrize("value", [1.0, -1.0, float("nan"), float("inf"), True])
@pytest.mark.parametrize("name", ["lambda_tasep_shock", "lambda_ot_cadence"])
def test_unused_objectives_rejected(name, value):
    with pytest.raises(ValueError):
        DPV240Optimizer(**{name: value})


def test_prefix_proxy_is_explicit():
    result = TransferMatrixFeasibilityGuard.check_feasibility(
        "KEL", AutomatonCompiler.compile(["GGTCTC"], include_rc=True)
    )
    assert result["scope"] == "PREFIX_REACHABILITY_ONLY"
    assert result["entropy_is_calibrated"] is False


@pytest.mark.parametrize(
    "protein,gc_min,gc_max",
    [("KK", 0.0, 1.0), ("KK", 0.15, 0.2), ("KEL", 0.2, 0.5), ("KEL", 0.0, 0.0)],
)
def test_dp_matches_independent_exhaustive_oracle(protein, gc_min, gc_max):
    """Enumerate every synonymous CDS independently of engine scoring/state helpers."""
    import math
    from itertools import product

    from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
    from factorforge.utils.exceptions import UnsatisfiableDesignError

    choices = [
        [c for c, aa in STANDARD_GENETIC_CODE.items() if aa == residue] for residue in protein
    ]
    # Artificial nonuniform weights, no host table or private sequence.
    weights = {c: 0.2 + (i + 1) / 100 for i, c in enumerate(sorted(STANDARD_GENETIC_CODE))}
    feasible = []
    for tokens in product(*choices):
        dna = "".join(tokens)
        gc = sum(b in "GC" for b in dna)
        if math.ceil(len(dna) * gc_min) <= gc <= math.floor(len(dna) * gc_max):
            feasible.append((sum(math.log(weights[c]) for c in tokens[1:]), dna))
    opt = DPV240Optimizer(
        ramp_length_codons=1,
        initiation_gc_min=0.0,
        initiation_gc_max=1.0,
        alpha_harmonization=0.0,
        beta_mfe_proxy=0.0,
        auto_annotate_pauses=False,
        forbidden_motifs=[],
        homopolymer_max_run=100,
    )
    if not feasible:
        with pytest.raises(UnsatisfiableDesignError):
            opt.optimize(protein, weights, target_gc_min=gc_min, target_gc_max=gc_max)
    else:
        result = opt.optimize(protein, weights, target_gc_min=gc_min, target_gc_max=gc_max)
        best = max(score for score, _ in feasible)
        assert result["score"] == pytest.approx(best, abs=1e-12)
        assert result["sequence"] in {dna for score, dna in feasible if abs(score - best) < 1e-12}
        assert (
            "".join(
                STANDARD_GENETIC_CODE[result["sequence"][i : i + 3]]
                for i in range(0, len(result["sequence"]), 3)
            )
            == protein
        )


def test_cadence_proxy_has_hand_calculated_value():
    from factorforge.engines.dp_v2_4 import WassersteinCadenceEvaluator

    # Candidate cumulative values: [2,3]; target: [1,3]; mean absolute difference=0.5.
    value = WassersteinCadenceEvaluator.evaluate_1d_wasserstein(
        "AAAAAG", {"AAA": 0.5, "AAG": 1.0}, {2: [("toy", 1.0, 0.5)]}
    )
    assert value == 0.5


def test_ramp_preserves_input_and_clips_upstream():
    from copy import deepcopy

    from factorforge.engines.dp_v2_4 import TASEPKineticsEngine

    original = {1: [("toy", 3.0, 0.3)], 3: [("other", 2.0, 0.4)]}
    snapshot = deepcopy(original)
    result = TASEPKineticsEngine.compute_smoothed_pause_map("KEL", original, smoothing_window=3)
    assert original == snapshot
    assert set(result) == {1, 2, 3}
    assert result[1] == original[1] and result[3] == original[3]


def test_prefix_reachability_does_not_prove_suffix_feasibility():
    from factorforge.utils.exceptions import UnsatisfiableDesignError

    # Prefix M is feasible; all AAA/AAG codons for the later K are forbidden.
    opt = DPV240Optimizer(
        forbidden_motifs=["AAA", "AAG"], ramp_length_codons=1, auto_annotate_pauses=False
    )
    prefix = TransferMatrixFeasibilityGuard.check_feasibility("MK", opt.automaton, ramp_length=1)
    assert prefix["feasible"] is True
    with pytest.raises(UnsatisfiableDesignError):
        opt.optimize(
            "MK", {"ATG": 1.0, "AAA": 0.8, "AAG": 0.4}, target_gc_min=0.0, target_gc_max=1.0
        )
