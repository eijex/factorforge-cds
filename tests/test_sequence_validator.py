"""Tests for sequence validation utilities."""

import pytest

from factorforge.utils.exceptions import SequenceValidationError
from factorforge.utils.sequence_validator import (
    detect_sequence_type,
    validate_cds_output,
    validate_and_normalize,
    validate_dna_sequence,
    validate_protein_sequence,
)


def test_detect_dna():
    assert detect_sequence_type("ATGGCC") == "dna"
    assert detect_sequence_type("ATGGCCAAATTT") == "dna"


def test_detect_protein():
    assert detect_sequence_type("MKKGEL") == "protein"
    assert detect_sequence_type("DEFHIL") == "protein"


def test_detect_ambiguous():
    assert detect_sequence_type("MA") == "ambiguous"
    assert detect_sequence_type("AT") == "ambiguous"


def test_validate_and_normalize_dna():
    seq, type_ = validate_and_normalize("atggcc", "auto")
    assert seq == "ATGGCC"
    assert type_ == "dna"


def test_validate_and_normalize_protein():
    seq, type_ = validate_and_normalize("MKKGEL", "auto")
    assert seq == "MKKGEL"
    assert type_ == "protein"


def test_validate_empty_sequence():
    with pytest.raises(SequenceValidationError, match="Empty sequence"):
        validate_and_normalize("")


def test_validate_invalid_characters():
    with pytest.raises(SequenceValidationError, match="Invalid characters"):
        validate_and_normalize("ATGXXX")


def test_validate_dna_sequence():
    assert validate_dna_sequence("ATGGCC") == "ATGGCC"

    with pytest.raises(SequenceValidationError, match="Invalid DNA characters"):
        validate_dna_sequence("ATGMKKGEL")


def test_validate_protein_sequence():
    assert validate_protein_sequence("MKKGEL") == "MKKGEL"

    with pytest.raises(SequenceValidationError):
        validate_protein_sequence("123ABC")


def test_validate_cds_output_passes_normal_cds():
    result = validate_cds_output("MAF", "ATGGCTTTC")

    assert result == {"passed": True, "errors": [], "aa_identity": 1.0}


def test_validate_cds_output_fails_internal_stop():
    result = validate_cds_output("MAF", "ATGTAATTC")

    assert result["passed"] is False
    assert "internal_stop_codon" in result["errors"]


def test_validate_cds_output_fails_aa_mismatch():
    result = validate_cds_output("MAF", "ATGGCTTAC")

    assert result["passed"] is False
    assert result["aa_identity"] == pytest.approx(2 / 3)
    assert any(error.startswith("aa_mismatch") for error in result["errors"])


def test_validate_cds_output_fails_length_error():
    result = validate_cds_output("MA", "ATGG")

    assert result["passed"] is False
    assert "length_not_divisible_by_3" in result["errors"]
