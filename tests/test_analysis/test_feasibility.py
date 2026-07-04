"""Tests for exact synonymous CAI/GC feasibility analysis."""

from __future__ import annotations

from factorforge.analysis.feasibility import analyze_feasibility
from factorforge.analysis.metrics import calculate_gc


WEIGHTS = {
    "ATG": 1.0,
    "GCT": 0.5,
    "GCC": 1.0,
    "GCA": 0.4,
    "GCG": 0.8,
    "AAA": 1.0,
    "AAG": 0.9,
}


def test_feasibility_short_protein_known_bounds() -> None:
    result = analyze_feasibility("MA", WEIGHTS, target_cai=0.9, target_gc_low=50, target_gc_high=70)

    assert result["minimum_possible_gc"] == 50.0
    assert result["maximum_possible_gc"] == 66.66666666666666
    assert result["target"]["feasible"] is True


def test_impossible_gc_target_flagged_impossible() -> None:
    result = analyze_feasibility(
        "MK", WEIGHTS, target_cai=0.1, target_gc_low=90, target_gc_high=100
    )

    assert result["target"]["feasible"] is False
    assert result["target"]["max_cai_under_gc"] is None


def test_feasible_target_flagged_feasible() -> None:
    result = analyze_feasibility("MK", WEIGHTS, target_cai=0.9, target_gc_low=16, target_gc_high=34)

    assert result["target"]["feasible"] is True
    assert result["target"]["best_candidate"]["dna_sequence"] == "ATGAAA"


def test_max_cai_without_gc_is_at_least_strict_gc_max() -> None:
    result = analyze_feasibility(
        "MAK", WEIGHTS, target_cai=0.1, target_gc_low=40, target_gc_high=50
    )
    strict = result["target"]["max_cai_under_gc"]

    assert strict is not None
    assert result["maximum_achievable_cai_without_gc"] >= strict


def test_local_gc_reporting_for_best_candidate() -> None:
    result = analyze_feasibility(
        "MAK", WEIGHTS, target_cai=0.1, target_gc_low=30, target_gc_high=70
    )
    candidate = result["target"]["best_candidate"]

    assert candidate is not None
    assert candidate["first_region_gc"]["first_30nt_gc"] == calculate_gc(candidate["dna_sequence"])
    assert "gc_window_min" in candidate
    assert "gc_window_max" in candidate


def test_cai_authority_unresolved_without_reference_id() -> None:
    result = analyze_feasibility("MA", WEIGHTS, target_cai=0.9, target_gc_low=50, target_gc_high=70)

    assert result["cai_authority"] == {
        "reference_id": None,
        "reference_role": "cai_evaluation",
        "reference_relationship": "unresolved",
    }
    assert result["target"]["best_candidate"]["cai_authority"]["reference_role"] == "cai_evaluation"


def test_cai_authority_same_as_generation_reference_with_reference_id() -> None:
    result = analyze_feasibility(
        "MA",
        WEIGHTS,
        target_cai=0.9,
        target_gc_low=50,
        target_gc_high=70,
        codon_reference_id="test_generation_ref_v1",
    )
    authority = result["target"]["best_candidate"]["cai_authority"]

    assert authority == {
        "reference_id": "test_generation_ref_v1",
        "reference_role": "cai_evaluation",
        "reference_relationship": "same_as_generation_reference",
    }
    assert result["cai_authority"]["reference_relationship"] != "distinct_from_generation_reference"


def test_dp_cai_authority_candidate_dicts_are_defensive_copies() -> None:
    result = analyze_feasibility(
        "MAK",
        WEIGHTS,
        target_cai=0.1,
        target_gc_low=30,
        target_gc_high=70,
        codon_reference_id="test_generation_ref_v1",
    )

    result["target"]["best_candidate"]["cai_authority"]["reference_id"] = "mutated"

    assert (
        result["best_candidate_without_gc"]["cai_authority"]["reference_id"]
        == "test_generation_ref_v1"
    )
