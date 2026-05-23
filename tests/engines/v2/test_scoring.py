"""
Unit tests for the multidimensional scoring module.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.scoring import (
    PROFILE_SCORING_CONFIGS,
    ScoringConfig,
    calculate_composite_score,
    normalize_mfe,
)


class TestScoringConfig:
    """Test ScoringConfig weight normalization."""

    def test_weights_normalize_to_one(self):
        """Active weights should sum to 1.0 after normalization."""
        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=True)
        total = cfg.w_cai + cfg.w_gc + cfg.w_mfe
        assert abs(total - 1.0) < 1e-6

    def test_mfe_disabled_renormalizes(self):
        """When use_mfe=False, w_mfe should be 0 and CAI+GC sum to 1."""
        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=False)
        assert cfg.w_mfe == 0.0
        total = cfg.w_cai + cfg.w_gc
        assert abs(total - 1.0) < 1e-6

    def test_profile_presets_exist(self):
        """All expected profile presets are available."""
        expected = {"balanced", "high_cai", "gc_target", "assembly_friendly", "ramp", "viral_delivery"}
        assert set(PROFILE_SCORING_CONFIGS.keys()) == expected


class TestCompositeScore:
    """Test calculate_composite_score function."""

    def test_perfect_cai_perfect_gc(self):
        """Perfect CAI and optimal GC should produce a high score."""
        score = calculate_composite_score(cai=1.0, gc=60.0, profile="balanced")
        assert score > 0.9

    def test_zero_cai(self):
        """Zero CAI should lower the score significantly."""
        score = calculate_composite_score(cai=0.0, gc=42.5, profile="balanced")
        assert score < 0.5

    def test_gc_far_from_optimum(self):
        """GC very far from optimum should lower the score."""
        score_good = calculate_composite_score(cai=0.8, gc=42.5, profile="balanced")
        score_bad = calculate_composite_score(cai=0.8, gc=80.0, profile="balanced")
        assert score_good > score_bad

    def test_high_cai_profile_weights(self):
        """High-CAI profile weights CAI more heavily."""
        # Use GC off-optimal so gc_target (w_gc=0.7) is penalized more than high_cai (w_gc=0.1)
        score_hcai = calculate_composite_score(cai=0.95, gc=42.5, profile="high_cai")
        score_gc = calculate_composite_score(cai=0.95, gc=42.5, profile="gc_target")
        # High-CAI should be higher because it weights CAI at 0.8 vs gc_target's 0.1
        assert score_hcai > score_gc

    def test_gc_target_profile_custom_gc(self):
        """GC-target profile uses custom target_gc."""
        score_at_target = calculate_composite_score(cai=0.5, gc=50.0, profile="gc_target", target_gc=50.0)
        score_off_target = calculate_composite_score(cai=0.5, gc=30.0, profile="gc_target", target_gc=50.0)
        assert score_at_target > score_off_target

    def test_unknown_profile_falls_back_to_balanced(self):
        """Unknown profile name falls back to balanced config."""
        score = calculate_composite_score(cai=0.8, gc=42.5, profile="nonexistent_profile")
        expected = calculate_composite_score(cai=0.8, gc=42.5, profile="balanced")
        assert score == expected

    def test_explicit_config_overrides_profile(self):
        """Explicit ScoringConfig overrides profile."""
        cfg = ScoringConfig(w_cai=1.0, w_gc=0.0, w_mfe=0.0, use_mfe=False)
        score = calculate_composite_score(cai=0.75, gc=0.0, config=cfg)
        assert abs(score - 0.75) < 0.01

    def test_score_clamped_0_to_1(self):
        """Score should always be between 0 and 1."""
        for cai in [0.0, 0.5, 1.0]:
            for gc in [0.0, 42.5, 100.0]:
                score = calculate_composite_score(cai=cai, gc=gc, profile="balanced")
                assert 0.0 <= score <= 1.0

    def test_ramp_profile_exists(self):
        """RAMP profile is recognized and produces valid score."""
        score = calculate_composite_score(cai=0.7, gc=42.5, profile="ramp")
        assert 0.0 <= score <= 1.0


class TestNormalizeMFE:
    """Test MFE normalization."""

    def test_zero_mfe_gives_one(self):
        """MFE = 0 (no structure) → normalized = 1.0."""
        assert normalize_mfe(0.0, 300) == 1.0

    def test_very_negative_mfe_gives_zero(self):
        """Very negative MFE → normalized ≈ 0.0."""
        result = normalize_mfe(-150.0, 300)  # -0.5 kcal/mol/nt
        assert abs(result) < 0.01

    def test_moderate_mfe(self):
        """Moderate MFE → intermediate value."""
        result = normalize_mfe(-75.0, 300)  # -0.25 kcal/mol/nt
        assert 0.4 < result < 0.6

    def test_zero_length_returns_neutral(self):
        """Zero-length sequence → neutral 0.5."""
        assert normalize_mfe(0.0, 0) == 0.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
