"""
Unit tests for RuleEngine
"""

import sys
from pathlib import Path

import pytest

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.rules.rule_engine import RuleEngine


@pytest.fixture
def sample_codon_table():
    """Minimal codon table for deterministic tests"""
    return {
        "amino_acids": {
            "M": {"codons": ["ATG"]},
            "N": {"codons": ["AAT", "AAC"]},
            "K": {"codons": ["AAA", "AAG"]},
            "A": {"codons": ["GCC", "GCT"]},
            "*": {"codons": ["TAA", "TAG", "TGA"]},
        },
        "codons": {
            "ATG": {"aa": "M", "frequency": 1.0},
            "AAT": {"aa": "N", "frequency": 0.5},
            "AAC": {"aa": "N", "frequency": 0.5},
            "AAA": {"aa": "K", "frequency": 0.5},
            "AAG": {"aa": "K", "frequency": 0.5},
            "GCC": {"aa": "A", "frequency": 0.5},
            "GCT": {"aa": "A", "frequency": 0.5},
            "TAA": {"aa": "*", "frequency": 0.34},
            "TAG": {"aa": "*", "frequency": 0.33},
            "TGA": {"aa": "*", "frequency": 0.33},
        },
    }


@pytest.fixture
def engine(sample_codon_table):
    """Create RuleEngine instance"""
    return RuleEngine(sample_codon_table)


class TestPolyAScanner:
    """Tests for PolyA scanning"""

    def test_scan_polya_single_pattern(self, engine):
        """Detect a single AATAAA signal"""
        sequence = "ATGAATAAAGCCTAA"
        violations = engine.scan_polya(sequence)

        assert any(v["type"] == "polya_signal" for v in violations)
        assert any(v["pattern"] == "AATAAA" for v in violations)
        assert any(v["position"] == 3 for v in violations)

    def test_scan_polya_multiple_patterns(self, engine):
        """Detect multiple PolyA patterns"""
        sequence = "AATAAAATTAAAGTAAA"
        violations = engine.scan_polya(sequence)

        signals = [v["pattern"] for v in violations if v["type"] == "polya_signal"]
        assert "AATAAA" in signals
        assert "ATTAAA" in signals
        assert "AGTAAA" in signals

    def test_scan_polya_no_patterns(self, engine):
        """Return no violations when patterns are absent"""
        sequence = "ATGGCCGCTGCTTAA"
        violations = engine.scan_polya(sequence)

        assert violations == []

    def test_scan_polya_empty_sequence(self, engine):
        """Handle empty sequence without errors"""
        violations = engine.scan_polya("")

        assert violations == []


class TestGCExtremes:
    """Tests for GC extremes scanning"""

    def test_scan_gc_extremes_high_gc(self, engine):
        """Detect high-GC windows"""
        sequence = "G" * 10
        violations = engine.scan_gc_extremes(sequence, window=5, min_gc=20, max_gc=80)

        assert len(violations) > 0
        assert all(v["type"] == "gc_extreme" for v in violations)
        assert all(v["gc"] == 100.0 for v in violations)


class TestFixViolation:
    """Tests for fix_violation"""

    def test_fix_violation_removes_polya(self, engine):
        """Fix a PolyA signal via synonymous substitution"""
        sequence = "ATGAATAAAGCCTAA"
        violations = engine.scan_polya(sequence)
        polya = next(v for v in violations if v["type"] == "polya_signal")

        result = engine.fix_violation(sequence, polya)

        assert result["success"] is True
        assert "AATAAA" not in result["modified_seq"]
        assert result["aa_preserved"] is True

    def test_fix_violation_no_synonyms(self):
        """Fail when no synonymous codons are available"""
        codon_table = {
            "amino_acids": {
                "M": {"codons": ["ATG"]},
                "N": {"codons": ["AAT"]},
                "K": {"codons": ["AAA"]},
                "*": {"codons": ["TAA"]},
            },
            "codons": {
                "ATG": {"aa": "M", "frequency": 1.0},
                "AAT": {"aa": "N", "frequency": 1.0},
                "AAA": {"aa": "K", "frequency": 1.0},
                "TAA": {"aa": "*", "frequency": 1.0},
            },
        }
        engine = RuleEngine(codon_table)
        sequence = "ATGAATAAATAA"
        violations = engine.scan_polya(sequence)
        polya = next(v for v in violations if v["type"] == "polya_signal")

        result = engine.fix_violation(sequence, polya)

        assert result["success"] is False
        assert result["modified_seq"] == sequence

    def test_fix_violation_invalid_length(self, engine):
        """Fail when sequence length is not divisible by 3"""
        sequence = "ATGAATAAA" + "G"  # length not divisible by 3
        violations = engine.scan_polya(sequence)
        polya = next(v for v in violations if v["type"] == "polya_signal")

        result = engine.fix_violation(sequence, polya)

        assert result["success"] is False
        assert result["aa_preserved"] is False


class TestScanModes:
    """Tests for scan_all mode/include/exclude behavior."""

    def test_scan_all_fast_mode_subset(self, engine):
        sequence = "ATGAATAAAGCCATTTA" + ("GCGCGC" * 5)
        full_results = engine.scan_all(sequence, mode="full")
        fast_results = engine.scan_all(sequence, mode="fast")

        assert "repeats" in full_results
        assert "dinucleotides" in full_results
        assert "repeats" not in fast_results
        assert "dinucleotides" not in fast_results
        assert "polya" in fast_results
        assert "gc_extremes" in fast_results

    def test_scan_all_with_include_exclude(self, engine):
        sequence = "ATGAATAAAGCCATTTA"
        results = engine.scan_all(
            sequence,
            mode="full",
            include=["polya", "are", "gc_extremes"],
            exclude=["are"],
        )
        assert set(results.keys()) == {"polya", "gc_extremes"}

    def test_scan_all_unknown_mode_raises(self, engine):
        with pytest.raises(ValueError, match="Unknown scan mode"):
            engine.scan_all("ATG", mode="turbo")

    def test_scan_all_unknown_scanner_raises(self, engine):
        with pytest.raises(ValueError, match="Unknown scanners"):
            engine.scan_all("ATG", include=["polya", "unknown_scanner"])
