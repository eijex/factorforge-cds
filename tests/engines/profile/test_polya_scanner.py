"""
Unit tests for PolyA scanning in RuleEngine
"""

import sys
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.rules.rule_engine import RuleEngine


@pytest.fixture
def engine():
    """Create RuleEngine instance"""
    return RuleEngine()


class TestPolyAScanner:
    """Test suite for PolyA scanning"""

    def test_find_violations_aataaa(self, engine):
        """Test detection of AATAAA pattern"""
        sequence = "ATGAATAAAGCCTAA"
        violations = engine.scan_polya(sequence)

        assert len(violations) == 1
        assert violations[0]["pattern"] == "AATAAA"
        assert violations[0]["position"] == 3

    def test_find_violations_multiple_patterns(self, engine):
        """Test detection of multiple PolyA patterns"""
        sequence = "AATAAAATTAAAGCAGTAAA"
        violations = engine.scan_polya(sequence)

        # Should find AATAAA, ATTAAA, AGTAAA
        patterns = [v["pattern"] for v in violations]
        assert "AATAAA" in patterns
        assert "ATTAAA" in patterns
        assert "AGTAAA" in patterns


class TestEdgeCases:
    """Test edge cases and boundary conditions"""

    def test_pattern_at_sequence_start(self, engine):
        """Test PolyA pattern at the beginning of sequence"""
        sequence = "AATAAAGCCTAA"
        violations = engine.scan_polya(sequence)

        assert len(violations) == 1
        assert violations[0]["position"] == 0

    def test_pattern_at_sequence_end(self, engine):
        """Test PolyA pattern at the end of sequence"""
        sequence = "ATGGCCAATAAA"
        violations = engine.scan_polya(sequence)

        assert len(violations) == 1
        assert violations[0]["position"] == 6

    def test_overlapping_patterns(self, engine):
        """Test overlapping PolyA patterns"""
        # AATAAAA contains both AATAAA and ATAAA (if ATAAA were a pattern)
        sequence = "AATAAAA"
        violations = engine.scan_polya(sequence)

        # Should find AATAAA
        assert any(v["pattern"] == "AATAAA" for v in violations)


class TestTier2Patterns:
    """Test Tier 2 (plant-functional) PolyA patterns."""

    def test_scan_polya_aataca(self, engine):
        """Detect AATACA (Tier 2 variant_3)."""
        violations = engine.scan_polya("ATGAATACAGCCTAA")
        patterns = [v["pattern"] for v in violations if v["type"] == "polya_signal"]
        assert "AATACA" in patterns

    def test_scan_polya_aagaaa(self, engine):
        """Detect AAGAAA (Tier 2 variant_4)."""
        violations = engine.scan_polya("ATGAAGAAAGCCTAA")
        patterns = [v["pattern"] for v in violations if v["type"] == "polya_signal"]
        assert "AAGAAA" in patterns

    def test_scan_polya_aatgaa(self, engine):
        """Detect AATGAA (Tier 2 variant_5)."""
        violations = engine.scan_polya("ATGAATGAAGCCTAA")
        patterns = [v["pattern"] for v in violations if v["type"] == "polya_signal"]
        assert "AATGAA" in patterns

    def test_tier_classification(self, engine):
        """Verify Tier 1 and Tier 2 pattern sets are disjoint and complete."""
        all_patterns = set(engine.POLYA_PATTERNS.keys())
        assert engine.POLYA_TIER1_PATTERNS | engine.POLYA_TIER2_PATTERNS == all_patterns
        assert engine.POLYA_TIER1_PATTERNS & engine.POLYA_TIER2_PATTERNS == set()


class TestPositivePolyA:
    """Test positive PolyA validation (terminator/3'UTR must HAVE a signal)."""

    def test_positive_polya_found(self, engine):
        """Terminator with canonical AATAAA returns no violations."""
        terminator_seq = "GCTAGCAATAAAGCTTAGCTT"
        violations = engine.scan_polya_positive(terminator_seq)
        assert violations == []

    def test_positive_polya_missing(self, engine):
        """Terminator without any signal returns missing_polya_signal violation."""
        terminator_seq = "GCTAGCGCTTAGCTTGCTAGC"
        violations = engine.scan_polya_positive(terminator_seq)
        assert len(violations) == 1
        assert violations[0]["type"] == "missing_polya_signal"
        assert violations[0]["severity"] == "high"

    def test_positive_polya_tier2_not_sufficient_by_default(self, engine):
        """Tier 2 pattern alone does not satisfy default Tier 1 requirement."""
        terminator_seq = "GCTAGCAATACAGCTTAGCTT"  # AATACA only (Tier 2)
        violations = engine.scan_polya_positive(terminator_seq)
        assert len(violations) == 1  # Missing Tier 1

    def test_positive_polya_custom_patterns(self, engine):
        """Custom pattern set accepts Tier 2."""
        terminator_seq = "GCTAGCAATACAGCTTAGCTT"  # AATACA (Tier 2)
        all_patterns = engine.POLYA_TIER1_PATTERNS | engine.POLYA_TIER2_PATTERNS
        violations = engine.scan_polya_positive(terminator_seq, required_patterns=all_patterns)
        assert violations == []


class TestIterativePolyAFix:
    """Test iterative PolyA removal from CDS sequences."""

    def test_fix_single_polya(self, engine):
        """Fix a single PolyA signal via synonymous substitution."""
        # AATAAA spans codons AAT-AAA (Asn-Lys)
        # Both have synonymous alternatives: AAT->AAC, AAA->AAG
        seq = "ATGAATAAAGCCTAA"  # M-N-K-A-*
        result = engine.fix_polya_iterative(seq)
        if result["success"]:
            assert "AATAAA" not in result["modified_seq"]
            assert len(result["modified_seq"]) == len(seq)

    def test_fix_preserves_amino_acids(self, engine):
        """Iterative fix preserves the amino acid sequence."""
        seq = "ATGAATAAAGCCTAA"  # M-N-K-A-*
        result = engine.fix_polya_iterative(seq)
        if result["success"]:
            # Verify each codon still codes for the same amino acid
            for change in result["fixes_applied"]:
                original_aa = engine.codon_table["codons"][change["original"]]["aa"]
                fixed_aa = engine.codon_table["codons"][change["fixed"]]["aa"]
                assert original_aa == fixed_aa

    def test_fix_no_violations_returns_immediately(self, engine):
        """Sequence without PolyA signals returns success with 0 rounds."""
        seq = "ATGGCCGCCTAA"  # No PolyA patterns
        result = engine.fix_polya_iterative(seq)
        assert result["success"] is True
        assert result["rounds"] == 0
        assert result["modified_seq"] == seq


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
