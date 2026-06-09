"""Explicit IUPAC and protein alphabet validation tests."""

import pytest

from factorforge.io.validation import SequenceValidationError, validate_sequence


def test_dna_strict_normalizes_whitespace_and_case() -> None:
    assert validate_sequence("atg aaa\nccc", "dna_strict") == "ATGAAACCC"


def test_dna_strict_rejects_ambiguity_code() -> None:
    with pytest.raises(SequenceValidationError):
        validate_sequence("ATGN", "dna_strict")


def test_dna_iupac_accepts_ambiguity_codes() -> None:
    assert validate_sequence("ACGTRYSWKMBDHVN", "dna_iupac") == "ACGTRYSWKMBDHVN"


def test_protein_strict_accepts_standard_amino_acids() -> None:
    assert validate_sequence("ACDEFGHIKLMNPQRSTVWY", "protein_strict")


def test_protein_strict_rejects_extended_amino_acids() -> None:
    with pytest.raises(SequenceValidationError):
        validate_sequence("MXX", "protein_strict")


def test_protein_extended_accepts_explicit_extensions() -> None:
    assert validate_sequence("MXBZUO*", "protein_extended") == "MXBZUO*"


def test_error_does_not_dump_full_raw_sequence() -> None:
    raw = "A" * 100 + "!"
    with pytest.raises(SequenceValidationError) as exc_info:
        validate_sequence(raw, "dna_strict")
    assert raw not in str(exc_info.value)
    assert "sha256-prefix=" in str(exc_info.value)
