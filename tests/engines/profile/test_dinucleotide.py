"""
Unit tests for CpG/TpA dinucleotide avoidance.

Tests cover:
- Utility functions (count, ratio)
- Rule engine scanner (scan_dinucleotides)
- Scoring integration (calculate_dinucleotide_score, composite score)
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.utils import calculate_dinucleotide_ratio, count_dinucleotides
from factorforge.engines.profile.rules.rule_engine import RuleEngine
from factorforge.engines.profile.scoring import (
    ScoringConfig,
    calculate_composite_score,
    calculate_dinucleotide_score,
)


class TestCountDinucleotides:
    """Tests for count_dinucleotides utility."""

    def test_count_cpg_basic(self):
        """Count CpG dinucleotides in simple sequence."""
        assert count_dinucleotides("ACGACG", "CG") == 2

    def test_count_tpa_basic(self):
        """Count TpA dinucleotides."""
        assert count_dinucleotides("ATAATA", "TA") == 2

    def test_count_empty_sequence(self):
        """Empty sequence returns 0."""
        assert count_dinucleotides("", "CG") == 0

    def test_count_single_base(self):
        """Single base returns 0."""
        assert count_dinucleotides("A", "CG") == 0

    def test_count_case_insensitive(self):
        """Case should not matter."""
        assert count_dinucleotides("acgacg", "CG") == 2
        assert count_dinucleotides("ACGACG", "cg") == 2

    def test_count_no_matches(self):
        """Sequence without target dinucleotide."""
        assert count_dinucleotides("AATTAATT", "CG") == 0


class TestDinucleotideRatio:
    """Tests for calculate_dinucleotide_ratio."""

    def test_ratio_short_sequence(self):
        """Very short sequence returns 0.0."""
        assert calculate_dinucleotide_ratio("A", "CG") == 0.0

    def test_ratio_no_target_bases(self):
        """Sequence without target bases returns 0.0."""
        assert calculate_dinucleotide_ratio("AATTAATT", "CG") == 0.0

    def test_ratio_returns_float(self):
        """Ratio is always a float."""
        result = calculate_dinucleotide_ratio("ACGTACGT", "CG")
        assert isinstance(result, float)

    def test_ratio_cpg_rich(self):
        """CpG-rich sequence has ratio > 1."""
        # Repeating CG directly creates high CpG density
        seq = "CGCGCGCGCGCGCGCG"
        ratio = calculate_dinucleotide_ratio(seq, "CG")
        assert ratio > 1.0

    def test_ratio_cpg_suppressed(self):
        """Sequence avoiding CpG has low ratio."""
        # C and G present but never adjacent as CG
        seq = "GCGCGCGCGCGCGCGC"  # GC but no CG
        # Actually "GC" repeated has CG at junction
        # Use a sequence that truly avoids CG
        seq = "CCCCGGGGCCCCGGGG"
        ratio = calculate_dinucleotide_ratio(seq, "CG")
        assert ratio < 1.0


class TestDinucleotideScanner:
    """Tests for RuleEngine.scan_dinucleotides."""

    @pytest.fixture
    def engine(self):
        """Create RuleEngine instance."""
        return RuleEngine()

    def test_high_cpg_detected(self, engine):
        """Sequence with many CpGs triggers violation."""
        seq = "ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG"
        violations = engine.scan_dinucleotides(seq, window=50)
        cpg_violations = [v for v in violations if v["dinucleotide"] == "CpG"]
        assert len(cpg_violations) > 0

    def test_clean_sequence_no_cpg_violations(self, engine):
        """AT-only sequence has no CpG violations."""
        seq = "AATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT"
        violations = engine.scan_dinucleotides(seq, window=50)
        cpg_violations = [v for v in violations if v["dinucleotide"] == "CpG"]
        assert len(cpg_violations) == 0

    def test_tpa_detected(self, engine):
        """High TpA density triggers violation."""
        seq = "TATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"
        violations = engine.scan_dinucleotides(seq, window=50)
        tpa_violations = [v for v in violations if v["dinucleotide"] == "TpA"]
        assert len(tpa_violations) > 0

    def test_scan_all_includes_dinucleotides(self, engine):
        """scan_all() includes dinucleotides key."""
        result = engine.scan_all("ATG" * 20)
        assert "dinucleotides" in result

    def test_violation_structure(self, engine):
        """Violation dicts have required keys."""
        seq = "ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG"
        violations = engine.scan_dinucleotides(seq, window=50)
        if violations:
            v = violations[0]
            assert "type" in v
            assert v["type"] == "dinucleotide_hotspot"
            assert "dinucleotide" in v
            assert "position" in v
            assert "density" in v
            assert "severity" in v

    def test_short_sequence_handled(self, engine):
        """Sequence shorter than window is handled gracefully."""
        violations = engine.scan_dinucleotides("ACGACG", window=50)
        # Should not raise, may or may not have violations
        assert isinstance(violations, list)


class TestDinucleotideScoring:
    """Tests for dinucleotide scoring integration."""

    def test_dinuc_score_range(self):
        """Score is always in [0, 1]."""
        score = calculate_dinucleotide_score("ATGATGATGATGATGATG")
        assert 0.0 <= score <= 1.0

    def test_dinuc_score_short_sequence(self):
        """Very short sequence scores 1.0 (too short to evaluate)."""
        assert calculate_dinucleotide_score("ATG") == 1.0

    def test_dinuc_score_cpg_rich_lower(self):
        """CpG-rich sequence scores lower than clean sequence."""
        clean = calculate_dinucleotide_score("AATTAATTAATTAATTAATTAATT")
        rich = calculate_dinucleotide_score("ACGACGACGACGACGACGACGACG")
        assert rich < clean

    def test_composite_score_dinuc_zero_backward_compat(self):
        """w_dinuc=0 produces same score as default profile."""
        score_default = calculate_composite_score(cai=0.8, gc=42.5, profile="balanced")
        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, w_dinuc=0.0)
        score_explicit = calculate_composite_score(cai=0.8, gc=42.5, config=cfg)
        assert abs(score_default - score_explicit) < 0.01

    def test_composite_score_dinuc_nonzero(self):
        """w_dinuc > 0 affects composite score and stays in range."""
        cfg = ScoringConfig(w_cai=0.4, w_gc=0.3, w_mfe=0.0, w_dinuc=0.3, use_mfe=False)
        score = calculate_composite_score(
            cai=0.8, gc=42.5, sequence="ATGATGATGATGATGATG", config=cfg
        )
        assert 0.0 <= score <= 1.0

    def test_scoring_config_normalizes_with_dinuc(self):
        """ScoringConfig normalizes all 4 weights to sum to 1.0."""
        cfg = ScoringConfig(w_cai=0.4, w_gc=0.2, w_mfe=0.2, w_dinuc=0.2)
        total = cfg.w_cai + cfg.w_gc + cfg.w_mfe + cfg.w_dinuc
        assert abs(total - 1.0) < 1e-6
