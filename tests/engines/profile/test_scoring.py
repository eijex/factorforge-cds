"""
Unit tests for the multidimensional scoring module.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.scoring import (
    GC_DECAY_WIDTH,
    GC_OPT_MAX,
    GC_OPT_MIN,
    PROFILE_SCORING_CONFIGS,
    ScoringConfig,
    calculate_composite_score,
    calculate_dinucleotide_score,
    gc_band_score,
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
        expected = {
            "balanced",
            "high_cai",
            "gc_target",
            "assembly_friendly",
            "ramp",
            "viral_delivery",
            "ml_enhanced",
        }
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
        """GC inside band scores higher than GC far outside band."""
        score_good = calculate_composite_score(cai=0.8, gc=60.0, profile="balanced")
        score_bad = calculate_composite_score(cai=0.8, gc=85.0, profile="balanced")
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


class TestGCBandScore:
    """Verify gc_band_score band function behaviour."""

    def test_inside_band_scores_one(self):
        assert gc_band_score(60.0, 55.0, 65.0) == 1.0

    def test_lower_boundary_scores_one(self):
        assert gc_band_score(55.0, 55.0, 65.0) == 1.0

    def test_upper_boundary_scores_one(self):
        assert gc_band_score(65.0, 55.0, 65.0) == 1.0

    def test_above_band_linear_decay(self):
        # 10 pp above gc_max=65, decay_width=20 → 1 - 10/20 = 0.5
        assert gc_band_score(75.0, 55.0, 65.0, decay_width=20.0) == pytest.approx(0.5)

    def test_below_band_linear_decay(self):
        # 10 pp below gc_min=55, decay_width=20 → 1 - 10/20 = 0.5
        assert gc_band_score(45.0, 55.0, 65.0, decay_width=20.0) == pytest.approx(0.5)

    def test_beyond_decay_width_scores_zero(self):
        # 20 pp above gc_max=65, decay_width=20 → 0.0
        assert gc_band_score(85.0, 55.0, 65.0, decay_width=20.0) == 0.0

    def test_far_outside_clamped_to_zero(self):
        assert gc_band_score(10.0, 55.0, 65.0) == 0.0

    def test_default_constants_used(self):
        """Default gc_min/gc_max/decay_width match module constants."""
        cfg = ScoringConfig()
        assert cfg.gc_min == GC_OPT_MIN
        assert cfg.gc_max == GC_OPT_MAX
        assert cfg.gc_decay_width == GC_DECAY_WIDTH

    def test_composite_score_uses_band_not_point(self):
        """GC inside band → higher gc contribution than same distance from gc_opt."""
        # GC=55 is at the lower boundary (score=1.0 with band, was 0.9 with old /50 formula)
        score_at_min = calculate_composite_score(cai=0.8, gc=55.0, profile="balanced")
        score_at_mid = calculate_composite_score(cai=0.8, gc=60.0, profile="balanced")
        # Both inside band → same GC component → scores equal
        assert score_at_min == score_at_mid

    def test_gc_target_uses_narrow_band_around_target(self):
        """gc_target profile: GC on-target scores higher than off-target."""
        on_target = calculate_composite_score(
            cai=0.8, gc=50.0, profile="gc_target", target_gc=50.0
        )
        off_target = calculate_composite_score(
            cai=0.8, gc=70.0, profile="gc_target", target_gc=50.0
        )
        assert on_target > off_target


class TestDinucleotideScore:
    """Verify CpG/TpA scoring weights."""

    def test_plant_default_cpg_inactive(self):
        """plant default: CpG inactive, so CpG-rich sequences are not penalized."""
        score = calculate_dinucleotide_score("CGCCGCCGCCGCCGCCGCCG", cpg_weight=0.0)
        assert score == pytest.approx(1.0)

    def test_tpa_rich_sequence_is_penalized(self):
        """TpA remains active in the plant default."""
        score = calculate_dinucleotide_score("TATATATATATATATATATATA", cpg_weight=0.0)
        assert score < 0.5

    def test_default_call_uses_plant_weights(self):
        """Default call keeps CpG inactive and TpA active."""
        sequence = "CGCCGCCGCCGCCGCCGCCG"
        assert calculate_dinucleotide_score(sequence) == calculate_dinucleotide_score(
            sequence, cpg_weight=0.0, tpa_weight=1.0
        )

    def test_mammalian_opt_in_penalizes_cpg(self):
        """Mammalian opt-in can penalize both CpG and TpA."""
        clean = calculate_dinucleotide_score(
            "CCCCGGGGCCCCGGGG", cpg_weight=1.0, tpa_weight=1.0
        )
        rich = calculate_dinucleotide_score(
            "CGCGCGCGCGCGCGCG", cpg_weight=1.0, tpa_weight=1.0
        )
        assert rich < clean


class TestAssemblyFriendlyProfile:
    """Verify assembly_friendly scoring is distinct from balanced."""

    def test_assembly_friendly_differs_from_balanced(self):
        """assembly_friendly must have different weights than balanced."""
        af = PROFILE_SCORING_CONFIGS["assembly_friendly"]
        bal = PROFILE_SCORING_CONFIGS["balanced"]
        assert af.w_cai != bal.w_cai, "w_cai must differ"
        assert af.w_gc != bal.w_gc, "w_gc must differ"

    def test_assembly_friendly_lower_cai_weight(self):
        """assembly_friendly reduces CAI pressure vs balanced."""
        af = PROFILE_SCORING_CONFIGS["assembly_friendly"]
        bal = PROFILE_SCORING_CONFIGS["balanced"]
        assert af.w_cai < bal.w_cai

    def test_assembly_friendly_higher_gc_weight(self):
        """assembly_friendly raises GC scoring weight vs balanced."""
        af = PROFILE_SCORING_CONFIGS["assembly_friendly"]
        bal = PROFILE_SCORING_CONFIGS["balanced"]
        assert af.w_gc > bal.w_gc

    def test_assembly_friendly_weights_sum_to_one(self):
        """assembly_friendly weights normalize correctly."""
        af = PROFILE_SCORING_CONFIGS["assembly_friendly"]
        total = af.w_cai + af.w_gc + af.w_mfe
        assert abs(total - 1.0) < 1e-6

    def test_balanced_unchanged(self):
        """Modifying assembly_friendly must not affect balanced weights."""
        bal = PROFILE_SCORING_CONFIGS["balanced"]
        assert abs(bal.w_cai - 0.5) < 1e-6
        assert abs(bal.w_gc - 0.3) < 1e-6


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
