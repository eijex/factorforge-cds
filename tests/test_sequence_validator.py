"""Tests for sequence validation utilities."""

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from factorforge.utils.exceptions import SequenceValidationError
from factorforge.utils.sequence_validator import (
    detect_sequence_type,
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
