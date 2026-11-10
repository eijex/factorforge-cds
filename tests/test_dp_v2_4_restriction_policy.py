"""Regression checks for the experimental full-construct motif policy."""

import math
from itertools import product

import pytest

from factorforge.constraints.type_iis import TYPE_IIS_MOTIFS, resolve_restriction_site_motifs
from factorforge.engines.dp_v2_4 import DPV240Optimizer
from factorforge.utils.exceptions import UnsatisfiableDesignError

BOUNDS = {"target_gc_min": 0, "target_gc_max": 1}
NOTI = "GCGGCCGC"


@pytest.mark.parametrize("enzyme", ["NotI", "noti", " NotI-HF "])
def test_noti_name_and_alias_resolve(enzyme):
    assert "NotI" not in TYPE_IIS_MOTIFS
    assert resolve_restriction_site_motifs([enzyme]) == [NOTI]


def test_unknown_name_fails_loudly():
    with pytest.raises(ValueError, match="Unknown restriction enzyme"):
        DPV240Optimizer(forbidden_enzymes=["Notl"])


@pytest.mark.parametrize(
    ("protein", "upstream", "stop", "downstream"),
    [
        ("M", NOTI, "TAA", ""),
        ("A", "GCGGCC", "TAA", ""),
        ("M", "", "", "CGGCCGC"),
        ("M", "", "TAA", NOTI),
    ],
)
def test_noti_is_checked_through_full_context(protein, upstream, stop, downstream):
    with pytest.raises(UnsatisfiableDesignError):
        DPV240Optimizer().optimize(
            protein,
            {},
            upstream_context=upstream,
            stop_codon=stop,
            downstream_context=downstream,
            **BOUNDS,
        )


def test_noti_cross_codon_and_literal_override():
    args = {"upstream_context": "GCG", "stop_codon": "", **BOUNDS}
    settings = {
        "ramp_length_codons": 1,
        "initiation_gc_min": 0,
        "initiation_gc_max": 1,
        "r_ramp_start": 1,
        "beta_mfe_proxy": 0,
    }
    default = DPV240Optimizer(**settings)
    override = DPV240Optimizer(forbidden_motifs=[], **settings)
    # A fixed prefix and an internal match exercise codon transitions independently.
    weights = {"GCG": 0.95, "GCC": 1.0, "GCT": 0.01, "GCA": 0.01}
    result = default.optimize("AA", weights, **args)
    assert NOTI not in "GCG" + result["sequence"]
    assert NOTI in result["motif_policy"]["resolved_enzyme_motifs"]
    result = override.optimize("AA", weights, **args)
    assert NOTI in "GCG" + result["sequence"]
    assert result["motif_policy"]["literal_motifs"] == []


def test_rc_flag_and_stop_junction_are_respected():
    motif = "TTACAT"
    DPV240Optimizer(forbidden_motifs=[motif], include_reverse_complement=False).optimize(
        "M", {}, **BOUNDS
    )
    with pytest.raises(UnsatisfiableDesignError):
        DPV240Optimizer(forbidden_motifs=[" ttacat "]).optimize("M", {}, **BOUNDS)


def test_strict_gc_empty_band_and_optional_fallback():
    optimizer = DPV240Optimizer()
    args = {"target_gc_min": 0, "target_gc_max": 0}
    with pytest.raises(UnsatisfiableDesignError):
        optimizer.optimize("M", {}, **args)
    result = optimizer.optimize("M", {}, allow_nearest_infeasible=True, **args)
    assert result["gc_feasible"] is False


@pytest.mark.parametrize(
    ("delta", "preferred"),
    [(0.0, "GCA"), (5e-13, "GCA"), (2e-12, "GCT")],
)
def test_real_arithmetic_oracle_and_near_tie_policy(delta, preferred):
    weights = {"ATG": 1.0, "GCA": 1.0 - delta, "GCT": 1.0, "GCC": 0.2, "GCG": 0.2}
    scores = {
        "".join(codons): math.log(weights[codons[1]])
        for codons in product(("ATG",), ("GCA", "GCC", "GCG", "GCT"), ("ATG",))
    }
    engine = DPV240Optimizer(
        forbidden_motifs=[], ramp_length_codons=1, initiation_gc_min=0, initiation_gc_max=1
    )
    result = engine.optimize("MAM", weights, **BOUNDS)
    sequence = result["sequence"]
    assert sequence == "ATG" + preferred + "ATG"
    assert 0 <= max(scores.values()) - scores[sequence] <= 1e-12


def test_gc_midpoint_breaks_exact_terminal_score_tie():
    engine = DPV240Optimizer(
        forbidden_motifs=[], ramp_length_codons=1, initiation_gc_min=0, initiation_gc_max=1
    )
    # ATG contributes one GC. GCA/GCT contribute two and GCC/GCG three.
    # With bounds [2, 4], midpoint 3 favors GCA/GCT, then lexical rank selects GCA.
    result = engine.optimize(
        "MA",
        {"ATG": 1, "GCA": 1, "GCT": 1, "GCC": 1, "GCG": 1},
        target_gc_min=2 / 6,
        target_gc_max=4 / 6,
    )
    assert result["sequence"] == "ATGGCA"
