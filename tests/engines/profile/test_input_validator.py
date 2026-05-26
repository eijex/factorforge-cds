"""
Unit tests for InputValidator
"""

import sys
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.validator import InputValidator, SequenceType


@pytest.fixture
def validator():
    """Create validator instance"""
    return InputValidator()


class TestSequenceTypeDetection:
    """Test sequence type auto-detection"""

    def test_detect_dna(self, validator):
        """Test DNA sequence detection"""
        assert validator.detect_sequence_type("ATGCGATCG") == SequenceType.DNA
        assert validator.detect_sequence_type("ATG GCC TAA") == SequenceType.DNA
        assert validator.detect_sequence_type("atgcgatcg") == SequenceType.DNA

    def test_detect_protein(self, validator):
        """Test protein sequence detection"""
        assert validator.detect_sequence_type("MAKLFG") == SequenceType.PROTEIN
        assert validator.detect_sequence_type("MAK LFG *") == SequenceType.PROTEIN
        assert validator.detect_sequence_type("maklfg") == SequenceType.PROTEIN

    def test_detect_fasta(self, validator):
        """Test FASTA format detection"""
        fasta = ">header\nATGCGATCG"
        assert validator.detect_sequence_type(fasta) == SequenceType.FASTA

    def test_detect_unknown(self, validator):
        """Test unknown sequence detection"""
        assert validator.detect_sequence_type("") == SequenceType.UNKNOWN
        assert validator.detect_sequence_type("12345") == SequenceType.UNKNOWN
        assert validator.detect_sequence_type("ATG123") == SequenceType.UNKNOWN


class TestDNAValidation:
    """Test DNA sequence validation"""

    def test_valid_dna(self, validator):
        """Test valid DNA sequence"""
        result = validator.validate("ATGGCCTAA")

        assert result["type"] == "dna"
        assert result["valid"] is True
        assert result["level"] == "valid"
        assert len(result["errors"]) == 0
        assert result["metadata"]["length"] == 9
        assert result["metadata"]["gc_content"] > 0

    def test_dna_with_spaces(self, validator):
        """Test DNA with whitespace"""
        result = validator.validate("ATG GCC TAA")

        assert result["type"] == "dna"
        assert result["processed_sequence"] == "ATGGCCTAA"

    def test_dna_frame_error(self, validator):
        """Test DNA with frame error"""
        result = validator.validate("ATGGCC")  # 6 bp, multiple of 3
        assert len(result["warnings"]) == 0

        result = validator.validate("ATGGCCA")  # 7 bp, not multiple of 3
        assert result["level"] == "warning"
        assert any(w["code"] == "FRAME_ERROR" for w in result["warnings"])

    def test_dna_frame_error_auto_fix(self, validator):
        """Test DNA frame error with auto-fix"""
        result = validator.validate("ATGGCCA", auto_fix=True)

        assert result["processed_sequence"] == "ATGGCC"
        assert any(w.get("auto_fixed") for w in result["warnings"])

    def test_dna_internal_stop(self, validator):
        """Test DNA with internal STOP codon"""
        result = validator.validate("ATGTAAGCC")  # TAA in middle

        assert any(w["code"] == "INTERNAL_STOP" for w in result["warnings"])

    def test_dna_extreme_gc(self, validator):
        """Test DNA with extreme GC content"""
        # Very low GC
        result = validator.validate("ATATATATATAT")
        assert any(w["code"] == "EXTREME_GC" for w in result["warnings"])

        # Very high GC
        result = validator.validate("GCGCGCGCGCGC")
        assert any(w["code"] == "EXTREME_GC" for w in result["warnings"])

    def test_dna_ambiguous_bases(self, validator):
        """Test DNA with ambiguous bases"""
        result = validator.validate("ATGNCCTAA")

        assert result["level"] == "warning"
        assert any(w["code"] == "AMBIGUOUS_DNA" for w in result["warnings"])
        assert "N" in result["metadata"]["ambiguous_bases"]

    def test_dna_invalid_chars(self, validator):
        """Test DNA with invalid characters"""
        result = validator.validate("ATG123")

        # Numbers cause the type to be detected as unknown
        assert result["type"] == "unknown"
        assert result["valid"] is False
        assert result["level"] == "error"


class TestProteinValidation:
    """Test protein sequence validation"""

    def test_valid_protein(self, validator):
        """Test valid protein sequence"""
        result = validator.validate("MAKLFG")

        assert result["type"] == "protein"
        assert result["valid"] is True
        assert result["level"] == "valid"
        assert result["metadata"]["length"] == 6

    def test_protein_with_stop(self, validator):
        """Test protein with STOP codon"""
        result = validator.validate("MAKLFG*")

        assert result["valid"] is True
        assert result["metadata"]["has_stop"] is True
        assert result["metadata"]["stop_positions"] == [7]

    def test_protein_internal_stop(self, validator):
        """Test protein with internal STOP"""
        result = validator.validate("MAK*LFG")

        assert result["level"] == "warning"
        assert any(w["code"] == "INTERNAL_STOP" for w in result["warnings"])

    def test_protein_ambiguous_aa(self, validator):
        """Test protein with ambiguous amino acids"""
        result = validator.validate("MAKXLFG")

        assert result["level"] == "warning"
        assert any(w["code"] == "AMBIGUOUS_AA" for w in result["warnings"])
        assert "X" in result["metadata"]["ambiguous_aa"]

    def test_protein_invalid_chars(self, validator):
        """Test protein with invalid characters"""
        result = validator.validate("MAK123")

        # Numbers cause the type to be detected as unknown
        assert result["type"] == "unknown"
        assert result["valid"] is False
        assert result["level"] == "error"


class TestFASTAValidation:
    """Test FASTA format validation"""

    def test_valid_fasta_dna(self, validator):
        """Test valid FASTA with DNA sequence"""
        fasta = """>test_sequence
ATGGCCTAA"""
        result = validator.validate(fasta)

        assert result["type"] == "fasta"
        assert result["fasta_header"] == "test_sequence"
        assert result["processed_sequence"] == "ATGGCCTAA"

    def test_valid_fasta_protein(self, validator):
        """Test valid FASTA with protein sequence"""
        fasta = """>test_protein
MAKLFG"""
        result = validator.validate(fasta)

        assert result["type"] == "fasta"
        assert result["fasta_header"] == "test_protein"

    def test_fasta_multiline(self, validator):
        """Test FASTA with multiline sequence"""
        fasta = """>multiline_seq
ATGGCC
TAAGCC
TAA"""
        result = validator.validate(fasta)

        assert result["processed_sequence"] == "ATGGCCTAAGCCTAA"

    def test_invalid_fasta_no_header(self, validator):
        """Test invalid FASTA without header"""
        fasta = "ATGGCCTAA"
        result = validator.validate(fasta)

        # Should be detected as DNA, not FASTA
        assert result["type"] == "dna"


class TestEdgeCases:
    """Test edge cases and boundary conditions"""

    def test_empty_sequence(self, validator):
        """Test empty sequence"""
        result = validator.validate("")

        assert result["type"] == "unknown"

    def test_whitespace_only(self, validator):
        """Test whitespace-only sequence"""
        result = validator.validate("   \n\t  ")

        assert result["type"] == "unknown"

    def test_lowercase_input(self, validator):
        """Test lowercase input"""
        result = validator.validate("atggcctaa")

        assert result["processed_sequence"] == "ATGGCCTAA"

    def test_mixed_case(self, validator):
        """Test mixed case input"""
        result = validator.validate("AtGgCcTaA")

        assert result["processed_sequence"] == "ATGGCCTAA"


class TestMetadata:
    """Test metadata generation"""

    def test_dna_metadata(self, validator):
        """Test DNA metadata"""
        result = validator.validate("ATGGCCTAA")

        assert "length" in result["metadata"]
        assert "gc_content" in result["metadata"]
        assert "has_stop" in result["metadata"]
        assert "stop_positions" in result["metadata"]

    def test_protein_metadata(self, validator):
        """Test protein metadata"""
        result = validator.validate("MAKLFG*")

        assert "length" in result["metadata"]
        assert "has_stop" in result["metadata"]
        assert "stop_positions" in result["metadata"]

    def test_gc_content_calculation(self, validator):
        """Test GC content calculation"""
        # 50% GC
        result = validator.validate("ATGC")
        assert result["metadata"]["gc_content"] == 50.0

        # 100% GC
        result = validator.validate("GCGCGC")
        assert result["metadata"]["gc_content"] == 100.0

        # 0% GC
        result = validator.validate("ATATAT")
        assert result["metadata"]["gc_content"] == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
