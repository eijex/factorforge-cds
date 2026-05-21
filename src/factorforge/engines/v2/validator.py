"""
Input Validator for FactorForge v2
Input validation and preprocessing module (P0-1)
"""

from __future__ import annotations

import re
from enum import Enum
from typing import Any


class SequenceType(Enum):
    """Sequence type enum"""

    PROTEIN = "protein"
    DNA = "dna"
    FASTA = "fasta"
    UNKNOWN = "unknown"


class ValidationLevel(Enum):
    """Validation level"""

    VALID = "valid"
    WARNING = "warning"
    ERROR = "error"


class InputValidator:
    """
    Input sequence validation and preprocessing

    Features:
    - Auto-detect AA/DNA/FASTA
    - Real-time input validation
    - Handle frame errors, stop codons, non-standard AAs
    """

    # Standard amino acids (20 + STOP)
    STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY*")

    # Ambiguous amino acids
    AMBIGUOUS_AA = {
        "B": "Asx (Asn or Asp)",
        "Z": "Glx (Gln or Glu)",
        "X": "Xaa (Unknown)",
        "J": "Xle (Leu or Ile)",
        "U": "Sec (Selenocysteine)",
        "O": "Pyl (Pyrrolysine)",
    }

    # DNA bases
    DNA_BASES = set("ATGC")
    AMBIGUOUS_DNA = set("NRYSWKMBDHV")  # IUPAC ambiguity codes

    def __init__(self) -> None:
        """Initialize"""
        self.warnings: list[dict[str, Any]] = []
        self.errors: list[dict[str, Any]] = []

    def detect_sequence_type(self, sequence: str) -> SequenceType:
        """
        Auto-detect sequence type

        Args:
            sequence: Input sequence (may include whitespace)

        Returns:
            SequenceType enum

        Raises:
            None.

        Examples:
            >>> validator = InputValidator()
            >>> validator.detect_sequence_type("ATGC").value
            'dna'
        """
        # Remove whitespace and newlines
        clean_seq = re.sub(r"\s+", "", sequence).upper()

        if not clean_seq:
            return SequenceType.UNKNOWN

        # Check FASTA format (starts with '>')
        if sequence.strip().startswith(">"):
            return SequenceType.FASTA

        # Analyze character set
        unique_chars = set(clean_seq)

        # DNA: only ATGC or ATGC + IUPAC codes
        if unique_chars <= (self.DNA_BASES | self.AMBIGUOUS_DNA):
            return SequenceType.DNA

        # Protein: amino acid characters
        if unique_chars <= (self.STANDARD_AA | set(self.AMBIGUOUS_AA.keys())):
            return SequenceType.PROTEIN

        return SequenceType.UNKNOWN

    def validate(self, sequence: str, auto_fix: bool = False) -> dict[str, Any]:
        """
        Validate input sequence

        Args:
            sequence: Input sequence
            auto_fix: Whether to auto-fix

        Returns:
            {
                "type": "protein|dna|fasta|unknown",
                "valid": True/False,
                "level": "valid|warning|error",
                "warnings": [...],
                "errors": [...],
                "metadata": {...},
                "processed_sequence": "..."  # Preprocessed sequence
            }

        Raises:
            None.

        Examples:
            >>> validator = InputValidator()
            >>> result = validator.validate("ATG GCC TAA")
            >>> result["type"]
            'dna'
        """
        self.warnings = []
        self.errors = []

        # 1. Detect sequence type
        seq_type = self.detect_sequence_type(sequence)

        # 2. Type-specific validation
        if seq_type == SequenceType.FASTA:
            return self._validate_fasta(sequence, auto_fix)
        elif seq_type == SequenceType.DNA:
            return self._validate_dna(sequence, auto_fix)
        elif seq_type == SequenceType.PROTEIN:
            return self._validate_protein(sequence, auto_fix)
        else:
            self.errors.append(
                {
                    "code": "UNKNOWN_TYPE",
                    "message": "Unable to detect sequence type.",
                    "suggestion": "Enter a DNA (ATGC) or Protein (20 AA) sequence.",
                }
            )
            return self._build_result(seq_type, sequence)

    def _validate_fasta(self, sequence: str, auto_fix: bool) -> dict[str, Any]:
        """Validate FASTA format"""
        lines = sequence.strip().split("\n")

        if not lines[0].startswith(">"):
            self.errors.append(
                {"code": "INVALID_FASTA", "message": "Missing FASTA header (must start with '>')."}
            )
            return self._build_result(SequenceType.FASTA, sequence)

        # Extract header
        header = lines[0][1:].strip()

        # Extract sequence lines
        seq_lines = [
            line.strip() for line in lines[1:] if line.strip() and not line.startswith(">")
        ]
        seq_content = "".join(seq_lines)

        # Recursively validate sequence content
        seq_result = self.validate(seq_content, auto_fix)

        # Add FASTA metadata
        seq_result["fasta_header"] = header
        seq_result["type"] = SequenceType.FASTA.value

        return seq_result

    def _validate_dna(self, sequence: str, auto_fix: bool) -> dict[str, Any]:
        """Validate DNA sequence"""
        clean_seq = re.sub(r"\s+", "", sequence).upper()

        # 1. Check invalid bases
        invalid_chars = set(clean_seq) - (self.DNA_BASES | self.AMBIGUOUS_DNA)
        if invalid_chars:
            self.errors.append(
                {
                    "code": "INVALID_DNA_CHARS",
                    "message": f"Invalid DNA bases: {', '.join(invalid_chars)}",
                    "suggestion": "Use only ATGC or IUPAC codes.",
                }
            )

        # 2. Warn on ambiguous bases
        ambiguous_chars = set(clean_seq) & self.AMBIGUOUS_DNA
        if ambiguous_chars:
            self.warnings.append(
                {
                    "code": "AMBIGUOUS_DNA",
                    "message": f"Ambiguous bases found: {', '.join(ambiguous_chars)}",
                    "suggestion": "Consider replacing with exact bases.",
                }
            )

        # 3. Frame check (multiple of 3)
        if len(clean_seq) % 3 != 0:
            remainder = len(clean_seq) % 3
            self.warnings.append(
                {
                    "code": "FRAME_ERROR",
                    "message": f"Sequence length is not a multiple of 3 ({len(clean_seq)} bp).",
                    "suggestion": f"Remove or add the last {remainder} bp.",
                    "auto_fix_option": {
                        "trim_end": len(clean_seq) - remainder,
                        "trim_start": len(clean_seq) - 3 + remainder,
                    },
                }
            )

            if auto_fix:
                clean_seq = clean_seq[: len(clean_seq) - remainder]
                self.warnings[-1]["auto_fixed"] = True

        # 4. Stop codon check
        stop_codons = {"TAA", "TAG", "TGA"}
        stop_positions = []

        for i in range(0, len(clean_seq) - 2, 3):
            codon = clean_seq[i : i + 3]
            if codon in stop_codons:
                codon_pos = i // 3 + 1  # 1-indexed
                stop_positions.append(codon_pos)

        # Warn on internal stop codons
        if len(stop_positions) > 1 or (
            len(stop_positions) == 1 and stop_positions[0] != len(clean_seq) // 3
        ):
            self.warnings.append(
                {
                    "code": "INTERNAL_STOP",
                    "message": f"Internal stop codon found: positions {stop_positions}",
                    "suggestion": "Verify whether this is intended.",
                }
            )

        # 5. GC content
        gc_count = clean_seq.count("G") + clean_seq.count("C")
        gc_content = (gc_count / len(clean_seq) * 100) if clean_seq else 0

        # Warn on extreme GC content
        if gc_content < 30 or gc_content > 70:
            self.warnings.append(
                {
                    "code": "EXTREME_GC",
                    "message": f"GC content is extreme: {gc_content:.1f}%",
                    "suggestion": "Recommended range is 30-70%.",
                }
            )

        metadata = {
            "length": len(clean_seq),
            "gc_content": round(gc_content, 2),
            "has_stop": len(stop_positions) > 0,
            "stop_positions": stop_positions,
            "ambiguous_bases": list(ambiguous_chars),
        }

        return self._build_result(SequenceType.DNA, clean_seq, metadata)

    def _validate_protein(self, sequence: str, auto_fix: bool) -> dict[str, Any]:
        """Validate protein sequence"""
        clean_seq = re.sub(r"\s+", "", sequence).upper()

        # 1. Check invalid amino acids
        invalid_chars = set(clean_seq) - (self.STANDARD_AA | set(self.AMBIGUOUS_AA.keys()))
        if invalid_chars:
            self.errors.append(
                {
                    "code": "INVALID_AA_CHARS",
                    "message": f"Invalid amino acids: {', '.join(invalid_chars)}",
                    "suggestion": "Use only standard 20 AAs or * (STOP).",
                }
            )

        # 2. Warn on ambiguous AAs
        ambiguous_found = {}
        for aa in self.AMBIGUOUS_AA:
            if aa in clean_seq:
                ambiguous_found[aa] = self.AMBIGUOUS_AA[aa]

        if ambiguous_found:
            self.warnings.append(
                {
                    "code": "AMBIGUOUS_AA",
                    "message": f"Ambiguous amino acids found: {ambiguous_found}",
                    "suggestion": "Consider replacing with specific amino acids.",
                }
            )

        # 3. Internal STOP check
        stop_positions = [i + 1 for i, aa in enumerate(clean_seq) if aa == "*"]

        if len(stop_positions) > 1 or (
            len(stop_positions) == 1 and stop_positions[0] != len(clean_seq)
        ):
            self.warnings.append(
                {
                    "code": "INTERNAL_STOP",
                    "message": f"Internal stop codon found: positions {stop_positions}",
                    "suggestion": "Verify whether this is intended.",
                }
            )

        metadata = {
            "length": len(clean_seq),
            "has_stop": len(stop_positions) > 0,
            "stop_positions": stop_positions,
            "ambiguous_aa": ambiguous_found,
        }

        return self._build_result(SequenceType.PROTEIN, clean_seq, metadata)

    def _build_result(
        self,
        seq_type: SequenceType,
        processed_seq: str,
        metadata: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Build validation result"""
        # Determine validation level
        if self.errors:
            level = ValidationLevel.ERROR
            valid = False
        elif self.warnings:
            level = ValidationLevel.WARNING
            valid = True
        else:
            level = ValidationLevel.VALID
            valid = True

        return {
            "type": seq_type.value,
            "valid": valid,
            "level": level.value,
            "warnings": self.warnings,
            "errors": self.errors,
            "metadata": metadata or {},
            "processed_sequence": processed_seq,
        }


# --- Usage example ---
if __name__ == "__main__":
    import json

    validator = InputValidator()

    # Test case 1: DNA sequence (valid)
    print("=== Test 1: Valid DNA ===")
    result1 = validator.validate("ATG GCC AAA TAA")
    print(json.dumps(result1, indent=2, ensure_ascii=False))

    # Test case 2: DNA sequence (frame error)
    print("\n=== Test 2: DNA with frame error ===")
    result2 = validator.validate("ATGGCCAA", auto_fix=True)
    print(json.dumps(result2, indent=2, ensure_ascii=False))

    # Test case 3: Protein sequence (ambiguous AA)
    print("\n=== Test 3: Protein with ambiguous AA ===")
    result3 = validator.validate("MAKXLF*")
    print(json.dumps(result3, indent=2, ensure_ascii=False))

    # Test case 4: FASTA format
    print("\n=== Test 4: FASTA format ===")
    fasta_input = """>GFP_test
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG
GTGGTGCCCATCCTGGTCGAGCTGGACGGCGAC
TAA"""
    result4 = validator.validate(fasta_input)
    print(json.dumps(result4, indent=2, ensure_ascii=False))
