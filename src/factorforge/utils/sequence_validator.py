"""Sequence validation and type detection utilities."""

from typing import Literal, Tuple

from .exceptions import SequenceValidationError

DNA_CHARS = set("ATGCN")
AMBIGUOUS_DNA_CHARS = set("ATGCM")
PROTEIN_ONLY_CHARS = set("DEFHIKLPQRSVWY")
VALID_CHARS = set("ACDEFGHIKLMNPQRSTVWY*")
MIN_DNA_LEN = 6


def _clean_sequence(seq: str) -> str:
    return "".join(seq.upper().split())


def detect_sequence_type(seq: str) -> Literal["dna", "protein", "ambiguous"]:
    """
    Detect if input sequence is DNA or protein.

    Args:
        seq: Input sequence string

    Returns:
        "dna": Valid DNA (only ATGCN)
        "protein": Valid protein (contains non-DNA amino acids)
        "ambiguous": Could be either (only contains A/T/G/C/M)

    Examples:
        >>> detect_sequence_type("ATGGCC")
        'dna'
        >>> detect_sequence_type("MKKGEL")
        'protein'
        >>> detect_sequence_type("MA")
        'ambiguous'
    """
    seq_upper = _clean_sequence(seq)
    if not seq_upper:
        return "ambiguous"

    seq_chars = set(seq_upper)

    # Protein-only letters present.
    if seq_chars & PROTEIN_ONLY_CHARS:
        return "protein"

    # DNA-only letters present.
    if seq_chars <= DNA_CHARS:
        if len(seq_upper) >= MIN_DNA_LEN and len(seq_upper) % 3 == 0:
            return "dna"
        return "ambiguous"

    # Ambiguous DNA (IUPAC M code).
    if seq_chars <= AMBIGUOUS_DNA_CHARS:
        return "ambiguous"

    # Default to protein for other amino-acid characters.
    return "protein"


def validate_and_normalize(
    seq: str,
    expected_type: Literal["dna", "protein", "auto"] = "auto",
) -> Tuple[str, Literal["dna", "protein"]]:
    """
    Validate and normalize sequence with type detection.

    Args:
        seq: Input sequence
        expected_type: Expected type or "auto" for auto-detection

    Returns:
        (normalized_sequence, detected_type)

    Raises:
        SequenceValidationError: If sequence is invalid or type mismatch

    Examples:
        >>> validate_and_normalize("atggcc", "auto")
        ('ATGGCC', 'dna')
        >>> validate_and_normalize("MKKGEL", "protein")
        ('MKKGEL', 'protein')
    """
    seq_clean = _clean_sequence(seq)

    if not seq_clean:
        raise SequenceValidationError("Empty sequence provided")

    invalid_chars = set(seq_clean) - VALID_CHARS
    if invalid_chars:
        raise SequenceValidationError(
            f"Invalid characters in sequence: {', '.join(sorted(invalid_chars))}"
        )

    detected = detect_sequence_type(seq_clean)

    if expected_type == "auto":
        if detected == "ambiguous":
            if len(seq_clean) >= MIN_DNA_LEN and len(seq_clean) % 3 == 0:
                return seq_clean, "dna"
            return seq_clean, "protein"
        return seq_clean, detected

    if expected_type != detected and detected != "ambiguous":
        raise SequenceValidationError(
            f"Expected {expected_type} sequence but detected {detected}. "
            f"Sequence: {seq_clean[:20]}{'...' if len(seq_clean) > 20 else ''}"
        )

    return seq_clean, expected_type


def validate_dna_sequence(seq: str) -> str:
    """
    Validate DNA sequence.

    Args:
        seq: DNA sequence

    Returns:
        Normalized DNA sequence

    Raises:
        SequenceValidationError: If not valid DNA
    """
    seq_clean = _clean_sequence(seq)
    invalid = set(seq_clean) - DNA_CHARS
    if invalid:
        raise SequenceValidationError(f"Invalid DNA characters: {', '.join(sorted(invalid))}")

    return seq_clean


def validate_protein_sequence(seq: str) -> str:
    """
    Validate protein sequence.

    Args:
        seq: Protein sequence

    Returns:
        Normalized protein sequence

    Raises:
        SequenceValidationError: If not valid protein
    """
    seq_clean = _clean_sequence(seq)
    aa_chars = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid = set(seq_clean) - aa_chars
    if invalid:
        raise SequenceValidationError(
            f"Invalid amino acid characters: {', '.join(sorted(invalid))}"
        )

    return seq_clean
