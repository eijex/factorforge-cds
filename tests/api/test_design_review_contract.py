from __future__ import annotations

from api.optimize import handler
from factorforge.analysis.metrics import translate_dna


def _handler() -> handler:
    return object.__new__(handler)


def test_300bp_dna_uses_cds_semantics_and_preserves_translation() -> None:
    original = "ATG" + ("GCT" * 97) + "GAA" + "TAA"
    result = _handler().optimize_sequence(
        original,
        "balanced",
        False,
        False,
        False,
        objective="feasibility_best",
        host_profile="nbenthamiana",
        return_candidates=True,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )
    optimized = result["optimized_sequence"]
    assert result["validation"]["input_type"] == "cds"
    assert result["original_length"] == 300
    assert result["optimized_length"] == 300
    assert len(optimized) == 300
    assert translate_dna(optimized).rstrip("*") == translate_dna(original).rstrip("*")


def test_api_review_decision_and_manual_acceptance_are_separate() -> None:
    original = "ATG" + ("GCT" * 20) + "GGTCTC" + ("GCT" * 25) + "TAA"
    criteria = {
        "overall_gc": {"mode": "ignored"},
        "local_gc": {"mode": "ignored"},
        "type_iis": {"mode": "required"},
        "repeats": {"mode": "ignored"},
        "homopolymers": {"mode": "ignored"},
        "forbidden_motifs": {"mode": "ignored"},
        "cai": {"mode": "ignored"},
    }
    result = _handler().optimize_sequence(
        original,
        "balanced",
        False,
        False,
        False,
        objective="feasibility_best",
        host_profile="nbenthamiana",
        return_candidates=True,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )
    result["optimized_sequence"] = original
    result["recommended_candidate"]["dna_sequence"] = original
    reviewed = _handler().attach_design_review(
        result,
        input_sequence=original,
        acceptance_criteria=criteria,
        reviewer_disposition={"disposition": "accept_with_exception", "reason": "Reviewed"},
    )

    assert reviewed["automated_decision"] == "FAIL"
    assert reviewed["reviewer_disposition"]["final_state"] == "MANUALLY_ACCEPTED"
    assert reviewed["reviewer_disposition"]["automated_decision"] == "FAIL"
    assert reviewed["acceptance_criteria_snapshot"]
