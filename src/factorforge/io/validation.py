"""Explicit DNA/protein alphabet validation without raw-sequence leakage."""

from __future__ import annotations

import hashlib
import re
from typing import Final

VALIDATION_ALPHABETS: Final[dict[str, frozenset[str]]] = {
    "dna_strict": frozenset("ACGT"),
    "dna_iupac": frozenset("ACGTRYSWKMBDHVN"),
    "protein_strict": frozenset("ACDEFGHIKLMNPQRSTVWY"),
    "protein_extended": frozenset("ACDEFGHIKLMNPQRSTVWYXBZUO*"),
}


class SequenceValidationError(ValueError):
    """Raised when a sequence does not satisfy its explicit alphabet contract."""


def _fingerprint(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:12]


def validate_sequence(sequence: str, mode: str = "dna_strict") -> str:
    """Normalize and validate a sequence for the requested alphabet mode.

    Whitespace is removed and letters are uppercased. Error messages include a
    short preview and fingerprint, never the complete input sequence.
    """
    if mode not in VALIDATION_ALPHABETS:
        choices = ", ".join(sorted(VALIDATION_ALPHABETS))
        raise ValueError(f"Unknown validation mode {mode!r}; expected one of: {choices}")
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    normalized = re.sub(r"\s+", "", sequence).upper()
    if not normalized:
        raise SequenceValidationError("Sequence is empty after whitespace normalization")

    invalid = sorted(set(normalized) - VALIDATION_ALPHABETS[mode])
    if invalid:
        preview = normalized[:20]
        raise SequenceValidationError(
            f"Invalid symbols for {mode}: {invalid}; preview={preview!r}; "
            f"sha256-prefix={_fingerprint(normalized)}"
        )
    return normalized
