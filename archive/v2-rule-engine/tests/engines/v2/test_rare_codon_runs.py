"""
Unit tests for rare codon run scanning in RuleEngine.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.rules.rule_engine import RuleEngine


@pytest.fixture
def rare_codon_table():
    """Codon table where AAA is rare relative to AAG."""
    return {
        "amino_acids": {
            "M": {"codons": ["ATG"]},
            "K": {"codons": ["AAA", "AAG"]},
            "*": {"codons": ["TAA", "TAG", "TGA"]},
        },
        "codons": {
            "ATG": {"aa": "M", "frequency": 1.0},
            "AAA": {"aa": "K", "frequency": 0.1},
            "AAG": {"aa": "K", "frequency": 1.0},
            "TAA": {"aa": "*", "frequency": 0.34},
            "TAG": {"aa": "*", "frequency": 0.33},
            "TGA": {"aa": "*", "frequency": 0.33},
        },
    }


@pytest.fixture
def engine(rare_codon_table):
    """Create RuleEngine instance with deterministic rare codon weights."""
    return RuleEngine(rare_codon_table)


class TestRareCodonRuns:
    """Tests for RuleEngine.scan_rare_codon_runs."""

    def test_detects_three_codon_run(self, engine):
        """Run length 3 is flagged by default."""
        violations = engine.scan_rare_codon_runs("ATG" + ("AAA" * 3) + "ATG")

        assert len(violations) == 1
        violation = violations[0]
        assert violation["type"] == "rare_codon_run"
        assert violation["position"] == 3
        assert violation["codon_index"] == 1
        assert violation["run_length"] == 3
        assert violation["codons"] == ["AAA", "AAA", "AAA"]

    def test_two_codon_run_is_below_default_threshold(self, engine):
        """Run length 2 is not flagged when min_run is 3."""
        violations = engine.scan_rare_codon_runs("ATG" + ("AAA" * 2) + "ATG")

        assert violations == []

    def test_run_length_five_is_high_severity(self, engine):
        """Run length 5 is high severity."""
        violations = engine.scan_rare_codon_runs("AAA" * 5)

        assert len(violations) == 1
        assert violations[0]["run_length"] == 5
        assert violations[0]["severity"] == "high"

    @pytest.mark.parametrize("run_length", [3, 4])
    def test_run_lengths_three_and_four_are_medium_severity(self, engine, run_length):
        """Run lengths 3 and 4 are medium severity."""
        violations = engine.scan_rare_codon_runs("AAA" * run_length)

        assert len(violations) == 1
        assert violations[0]["run_length"] == run_length
        assert violations[0]["severity"] == "medium"

    def test_stop_codon_breaks_rare_codon_run(self, engine):
        """Stop codons are excluded and break a rare codon run."""
        violations = engine.scan_rare_codon_runs(("AAA" * 2) + "TAA" + ("AAA" * 2))

        assert violations == []

    def test_empty_and_non_codon_length_sequences_return_empty(self, engine):
        """Empty sequence and non-triplet sequence lengths return no violations."""
        assert engine.scan_rare_codon_runs("") == []
        assert engine.scan_rare_codon_runs("AAAA") == []

    def test_scan_all_full_mode_includes_rare_codon_runs(self, engine):
        """scan_all() includes rare_codon_runs in full mode."""
        result = engine.scan_all("ATG" * 5)

        assert "rare_codon_runs" in result

    def test_scan_all_fast_mode_excludes_rare_codon_runs(self, engine):
        """scan_all() excludes rare_codon_runs in fast mode."""
        result = engine.scan_all("ATG" * 5, mode="fast")

        assert "rare_codon_runs" not in result
