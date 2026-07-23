from __future__ import annotations

import pytest

from factorforge.analysis.metrics import translate_dna
from factorforge.design_review import (
    apply_reviewer_disposition,
    assert_pathway_invariants,
    default_acceptance_criteria,
    evaluate_candidate,
    normalize_acceptance_criteria,
    parse_sequence_input,
    restore_cds_stop_policy,
    original_cai,
)
from factorforge.analysis.metrics import load_codon_usage_table


BALANCED_DNA = "ATG" + ("GCT" * 97) + "GAA" + "TAA"


def test_plain_and_single_fasta_dna_are_normalized() -> None:
    plain = parse_sequence_input(BALANCED_DNA.lower())
    fasta = parse_sequence_input(f">fixture\n{BALANCED_DNA[:90]}\n{BALANCED_DNA[90:]}\n")

    assert plain["normalized_sequence"] == BALANCED_DNA
    assert fasta["normalized_sequence"] == BALANCED_DNA
    assert fasta["fasta_header"] == "fixture"
    assert fasta["input_type"] == "cds"


def test_multifasta_is_rejected_without_concatenation() -> None:
    with pytest.raises(ValueError, match="Multiple FASTA records detected"):
        parse_sequence_input(">one\nATGGCTTAA\n>two\nATGGAATAA")


def test_invalid_characters_and_frame_are_reported() -> None:
    with pytest.raises(ValueError, match="Invalid sequence characters: B"):
        parse_sequence_input("ATGBCC")

    parsed = parse_sequence_input("ATGG")
    assert parsed["input_type"] == "cds"
    assert parsed["generation_allowed"] is False
    assert parsed["summary"]["frame_valid"] is False


def test_protein_is_classified_separately() -> None:
    parsed = parse_sequence_input("MSKGEELFTGVV")
    assert parsed["input_type"] == "protein"
    assert parsed["summary"]["protein_length_aa"] == 12
    assert "reverse translate" in parsed["message"]


def test_300bp_cds_context_and_stop_policy_preserve_units() -> None:
    parsed = parse_sequence_input(BALANCED_DNA)
    candidate = restore_cds_stop_policy("ATG" + ("GCC" * 97) + "GAG", parsed)

    assert len(BALANCED_DNA) == 300
    assert parsed["summary"]["cds_length_bp"] == 300
    assert parsed["summary"]["codon_count"] == 100
    assert parsed["summary"]["protein_length_aa"] == 99
    assert len(candidate) == 300
    assert translate_dna(candidate).rstrip("*") == translate_dna(BALANCED_DNA).rstrip("*")
    assert_pathway_invariants(BALANCED_DNA, candidate, parsed)


def test_internal_stop_disables_generation() -> None:
    parsed = parse_sequence_input("ATGGCTTAAGCTGAA")
    assert parsed["generation_allowed"] is False
    assert parsed["summary"]["internal_stop_count"] == 1


def test_required_preferred_and_ignored_decisions() -> None:
    sequence = "ATG" + ("GCT" * 30) + "GGTCTC" + ("GCT" * 10) + "TAA"
    criteria = default_acceptance_criteria(0, 100)
    required = evaluate_candidate(sequence, cai=0.9, criteria=criteria)
    assert required["automated_decision"] == "FAIL"

    criteria["type_iis"]["mode"] = "preferred"
    preferred = evaluate_candidate(sequence, cai=0.9, criteria=criteria)
    assert preferred["automated_decision"] == "CONDITIONAL_PASS"

    criteria["type_iis"]["mode"] = "ignored"
    criteria["local_gc"]["mode"] = "ignored"
    criteria["repeats"]["mode"] = "ignored"
    criteria["homopolymers"]["mode"] = "ignored"
    criteria["forbidden_motifs"]["mode"] = "ignored"
    ignored = evaluate_candidate(sequence, cai=0.9, criteria=criteria)
    assert ignored["automated_decision"] == "PASS"


def test_preferred_repeat_issue_is_conditional() -> None:
    repeat = "ATGGCTGAACTGTTTGGT" * 2
    criteria = default_acceptance_criteria(0, 100)
    criteria["type_iis"]["mode"] = "ignored"
    criteria["local_gc"]["mode"] = "ignored"
    criteria["homopolymers"]["mode"] = "ignored"
    criteria["forbidden_motifs"]["mode"] = "ignored"
    result = evaluate_candidate(repeat, cai=0.9, criteria=criteria)
    assert result["automated_decision"] == "CONDITIONAL_PASS"
    assert result["preferred_warning_count"] >= 1


def test_manual_acceptance_of_fail_requires_reason_and_retains_fail() -> None:
    with pytest.raises(ValueError, match="written reason"):
        apply_reviewer_disposition("FAIL", {"disposition": "accept_with_exception"})

    result = apply_reviewer_disposition(
        "FAIL",
        {"disposition": "accept_with_exception", "reason": "Documented synthesis exception"},
    )
    assert result["final_state"] == "MANUALLY_ACCEPTED"
    assert result["automated_decision"] == "FAIL"


def test_criteria_validation_rejects_bad_mode() -> None:
    with pytest.raises(ValueError, match="must be required"):
        normalize_acceptance_criteria(
            {"cai": {"mode": "blocking"}},
            gc_min=40,
            gc_max=55,
        )


def test_required_type_iis_sites_and_forbidden_fixture_are_visible() -> None:
    fixture = "ATG" + "GGTCTC" + "CGTCTC" + "GCTCTTC" + "AATAAA" + ("A" * 10) + "TAA"
    criteria = default_acceptance_criteria(0, 100)
    result = evaluate_candidate(fixture, cai=0.9, criteria=criteria)
    assert result["automated_decision"] == "FAIL"
    assert len(result["details"]["type_iis_sites"]) >= 3
    assert result["details"]["forbidden_motifs"]
    assert result["details"]["longest_homopolymer"] >= 10


def test_original_cai_is_calculated_from_the_active_reference() -> None:
    table = load_codon_usage_table()
    high = "ATG" + "GCC" * 4 + "TAA"
    low = "ATG" + "GCT" * 4 + "TAA"
    high_cai = original_cai(high, table.codon_weights)
    low_cai = original_cai(low, table.codon_weights)
    assert high_cai is not None and low_cai is not None
    assert high_cai != low_cai


def test_all_criteria_pass_returns_pass() -> None:
    criteria = default_acceptance_criteria(0, 100)
    for criterion in criteria.values():
        criterion["mode"] = "ignored"
    result = evaluate_candidate("ATGGCTGAA", cai=0.9, criteria=criteria)
    assert result["automated_decision"] == "PASS"
