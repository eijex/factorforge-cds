"""Custom exceptions for FactorForge."""


class FactorForgeError(Exception):
    """Base exception for FactorForge."""


class SequenceValidationError(FactorForgeError):
    """Raised when sequence validation fails."""


class OptimizationError(FactorForgeError):
    """Raised when optimization fails."""


class EmptyCandidateError(OptimizationError):
    """Raised when no valid codon candidates are generated."""

    def __init__(self, amino_acid: str, reason: str = "") -> None:
        self.amino_acid = amino_acid
        message = f"No valid codon candidates for amino acid '{amino_acid}'"
        if reason:
            message = f"{message}: {reason}"
        super().__init__(message)


class FileFormatError(FactorForgeError):
    """Raised when file format is invalid."""


class CodonTableError(FactorForgeError):
    """Raised when codon table is invalid or missing."""
