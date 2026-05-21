"""Tests for structured candidate sequence validation."""

from __future__ import annotations

from factorforge.utils.validation import validate_candidate_sequence


def test_validator_passes_valid_sequence() -> None:
    result = validate_candidate_sequence("MAF", "ATGGCTTTC")

    assert result["passed"] is True
    assert result["amino_acid_identity"] == 1.0
    assert result["invalid_codon_count"] == 0
    assert result["internal_stop_count"] == 0


def test_validator_fails_sequence_not_divisible_by_three() -> None:
    result = validate_candidate_sequence("MA", "ATGG")

    assert result["passed"] is False
    assert "sequence length is not divisible by 3" in result["errors"]
    assert result["invalid_codon_count"] == 1


def test_validator_fails_internal_stop() -> None:
    result = validate_candidate_sequence("MAF", "ATGTAATTC")

    assert result["passed"] is False
    assert result["internal_stop_count"] == 1
    assert "internal stop codons detected" in result["errors"]


def test_validator_forbidden_motif_hard_fail() -> None:
    result = validate_candidate_sequence(
        "MGL",
        "ATGGGTCTC",
        config={"forbidden_motifs": ["GGTCTC"], "fail_forbidden_motifs": True},
    )

    assert result["passed"] is False
    assert result["forbidden_motif_count"] == 1
    assert "forbidden motifs detected" in result["errors"]


def test_validator_reports_homopolymers_repeats_and_gc_windows() -> None:
    result = validate_candidate_sequence(
        "MKK",
        "ATGAAAAAAAA",
        config={"gc_window_size": 6, "gc_window_step": 3, "gc_window_low": 40, "gc_window_high": 60},
    )

    assert result["passed"] is False
    assert result["homopolymer_count"] >= 1
    assert result["gc_window_outlier_count"] >= 1
    assert "local GC window outliers detected" in result["warnings"]

