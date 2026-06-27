"""
Unit tests for the multidimensional scoring module.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile import scoring
from factorforge.engines.profile.scoring import (
    GC_DECAY_WIDTH,
    GC_OPT_MAX,
    GC_OPT_MID,
    GC_OPT_MIN,
    MFE_MAX_SEQUENCE_LENGTH,
    PROFILE_SCORING_CONFIGS,
    ScoringConfig,
    calculate_composite_score,
    calculate_dinucleotide_score,
    calculate_mfe,
    compute_mfe_evidence,
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
        score = calculate_composite_score(cai=1.0, gc=GC_OPT_MID, profile="balanced")
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
        # Use GC off-optimal (outside the band, decay_width past GC_OPT_MAX) so
        # gc_target (w_gc=0.7) is penalized more than high_cai (w_gc=0.1).
        off_optimal_gc = GC_OPT_MAX + GC_DECAY_WIDTH / 2
        score_hcai = calculate_composite_score(cai=0.95, gc=off_optimal_gc, profile="high_cai")
        score_gc = calculate_composite_score(cai=0.95, gc=off_optimal_gc, profile="gc_target")
        # High-CAI should be higher because it weights CAI at 0.8 vs gc_target's 0.1
        assert score_hcai > score_gc

    def test_gc_target_profile_custom_gc(self):
        """GC-target profile uses custom target_gc."""
        score_at_target = calculate_composite_score(cai=0.5, gc=50.0, profile="gc_target", target_gc=50.0)
        score_off_target = calculate_composite_score(cai=0.5, gc=30.0, profile="gc_target", target_gc=50.0)
        assert score_at_target > score_off_target

    def test_unknown_profile_raises(self):
        """Unknown profile name raises, rather than silently using balanced."""
        with pytest.raises(ValueError, match="Unknown profile"):
            calculate_composite_score(cai=0.8, gc=42.5, profile="nonexistent_profile")

    def test_target_gc_rejected_for_non_gc_profile(self):
        """target_gc kwarg is only valid for the gc_target profile."""
        with pytest.raises(ValueError, match="target_gc"):
            calculate_composite_score(cai=0.8, gc=60.0, profile="balanced", target_gc=30.0)

    def test_negative_weight_rejected(self):
        """ScoringConfig rejects negative weights."""
        with pytest.raises(ValueError, match="w_gc"):
            ScoringConfig(w_cai=2.0, w_gc=-1.0, w_mfe=0.0, use_mfe=False)

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
        """MFE = 0 (no structure) → above the [-0.40, -0.15] clamp ceiling → normalized = 1.0."""
        assert normalize_mfe(0.0, 300) == 1.0

    def test_very_negative_mfe_gives_zero(self):
        """Very negative MFE → below the [-0.40, -0.15] clamp floor → normalized ≈ 0.0."""
        result = normalize_mfe(-150.0, 300)  # -0.5 kcal/mol/nt, below the -0.40 floor
        assert abs(result) < 0.01

    def test_moderate_mfe(self):
        """Moderate MFE at the clamp midpoint → intermediate value."""
        result = normalize_mfe(-82.5, 300)  # -0.275 kcal/mol/nt, midpoint of [-0.40, -0.15]
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
        # GC_OPT_MIN is at the lower boundary; GC_OPT_MID is the band midpoint.
        # Both fall inside [GC_OPT_MIN, GC_OPT_MAX] (score=1.0 with band scoring).
        score_at_min = calculate_composite_score(cai=0.8, gc=GC_OPT_MIN, profile="balanced")
        score_at_mid = calculate_composite_score(cai=0.8, gc=GC_OPT_MID, profile="balanced")
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


class TestMFEViennaBranches:
    """Cover the ViennaRNA-present / ViennaRNA-absent branches that
    normalize_mfe-only tests skip (most other tests call with
    sequence=None, bypassing the MFE branch entirely)."""

    SEQ = "ATG" + "GCT" * 30 + "TAA"

    def test_mfe_computed_contributes_to_score(self, monkeypatch):
        """With ViennaRNA available, MFE weight actually changes the score
        relative to an otherwise-identical MFE-disabled config."""
        monkeypatch.setattr(scoring, "_check_vienna_available", lambda: True)
        monkeypatch.setattr(scoring, "calculate_mfe", lambda seq: -40.0)

        cfg_with_mfe = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=True)
        cfg_without_mfe = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=False)

        score_with = calculate_composite_score(
            cai=0.5, gc=60.0, sequence=self.SEQ, config=cfg_with_mfe
        )
        score_without = calculate_composite_score(
            cai=0.5, gc=60.0, sequence=self.SEQ, config=cfg_without_mfe
        )
        assert score_with != score_without

    def test_vienna_unavailable_sets_not_computed(self, monkeypatch):
        """When ViennaRNA is unavailable, evidence reports not_computed
        honestly rather than a misleading 0.0."""
        monkeypatch.setattr(scoring, "_check_vienna_available", lambda: False)
        ev = compute_mfe_evidence(self.SEQ, profile="balanced")
        assert ev["mfe_kcal_mol"] is None
        assert ev["mfe_status"] == "not_computed"
        assert ev["mfe_used"] is False
        assert ev["mfe_warning"] is not None

    def test_mfe_failure_removes_weight(self, monkeypatch):
        """If calculate_mfe raises/returns None, the MFE weight is dropped
        (renormalized) rather than treated as a real 0.0 kcal/mol value."""
        monkeypatch.setattr(scoring, "_check_vienna_available", lambda: True)
        monkeypatch.setattr(scoring, "calculate_mfe", lambda seq: None)

        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=True)
        score = calculate_composite_score(cai=1.0, gc=60.0, sequence=self.SEQ, config=cfg)
        # With MFE dropped, score should equal CAI+GC-only scoring (renormalized).
        cfg_no_mfe = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=False)
        score_no_mfe = calculate_composite_score(
            cai=1.0, gc=60.0, sequence=self.SEQ, config=cfg_no_mfe
        )
        assert score == pytest.approx(score_no_mfe)

    def test_effective_weights_are_renormalized(self):
        """CAI/GC weights scale up to compensate when MFE is disabled,
        and still sum to 1.0."""
        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=False)
        assert cfg.w_mfe == 0.0
        assert abs((cfg.w_cai + cfg.w_gc) - 1.0) < 1e-6

    def test_use_mfe_false_never_calls_vienna(self, monkeypatch):
        """use_mfe=False must never invoke ViennaRNA at all, not just
        zero-weight it after computing."""
        called = []
        monkeypatch.setattr(scoring, "_check_vienna_available", lambda: True)
        monkeypatch.setattr(
            scoring, "calculate_mfe", lambda seq: called.append(seq) or -40.0
        )

        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.0, use_mfe=False)
        calculate_composite_score(cai=0.5, gc=60.0, sequence=self.SEQ, config=cfg)
        assert called == []


class TestMFELengthGuard:
    """170-fix: calculate_mfe() must not attempt RNA.fold() on sequences long
    enough to make ViennaRNA's O(n^3) Zuker MFE algorithm a multi-second-to-
    multi-minute cost (confirmed via faulthandler + an isolated RNA.fold()
    timing curve: 1000nt~2.2s, 2000nt~9.8s, 3000nt~24.6s on the dev machine) —
    reachable by an unauthenticated request within the existing public API
    length limits (5000aa/15000bp), i.e. an algorithmic-complexity DoS
    (CWE-407) before this guard existed."""

    def test_sequence_at_cutoff_is_not_skipped(self, monkeypatch):
        """Exactly MFE_MAX_SEQUENCE_LENGTH must still attempt the fold."""
        monkeypatch.setattr(scoring, "_check_vienna_available", lambda: True)
        import types
        monkeypatch.setitem(
            sys.modules, "RNA", types.SimpleNamespace(fold=lambda seq: (None, -10.0))
        )

        seq = "ATG" * (MFE_MAX_SEQUENCE_LENGTH // 3)
        result = calculate_mfe(seq[:MFE_MAX_SEQUENCE_LENGTH])
        assert result == -10.0

    def test_sequence_over_cutoff_returns_none_without_folding(self, monkeypatch, caplog):
        """One nucleotide over the cutoff must skip RNA.fold() entirely and
        log a warning explaining why (so a future 'MFE score looks wrong for
        long sequences' bug report has an immediate, discoverable answer)."""
        monkeypatch.setattr(scoring, "_check_vienna_available", lambda: True)

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError("RNA.fold() must not be called above the length cutoff")

        import types
        monkeypatch.setitem(
            sys.modules, "RNA", types.SimpleNamespace(fold=fail_if_called)
        )

        seq = "A" * (MFE_MAX_SEQUENCE_LENGTH + 1)
        with caplog.at_level("WARNING"):
            result = calculate_mfe(seq)
        assert result is None
        assert any(
            "MFE_MAX_SEQUENCE_LENGTH" in record.message for record in caplog.records
        )

    def test_compute_mfe_evidence_length_skip_is_not_mislabeled_as_failure(self):
        """170-fix follow-up (Opus review): a length-based skip must not be
        reported via the generic 'MFE computation failed' message — that
        wording implies a fold error, which would misdirect a caller
        debugging why MFE is missing for a long sequence."""
        seq = "A" * (MFE_MAX_SEQUENCE_LENGTH + 1)
        ev = compute_mfe_evidence(seq, profile="balanced")
        assert ev["mfe_used"] is False
        assert ev["mfe_warning"] is not None
        assert "failed" not in ev["mfe_warning"].lower()
        assert str(MFE_MAX_SEQUENCE_LENGTH) in ev["mfe_warning"]

    def test_over_cutoff_drops_mfe_weight_like_vienna_unavailable(self):
        """Above the cutoff, calculate_composite_score must behave exactly
        like the pre-existing 'ViennaRNA unavailable' path (weight dropped,
        renormalized) — no new state, just reusing the existing fallback."""
        long_seq = "ATG" + "GCT" * ((MFE_MAX_SEQUENCE_LENGTH // 3) + 1) + "TAA"
        cfg = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=True)
        score = calculate_composite_score(cai=1.0, gc=60.0, sequence=long_seq, config=cfg)

        cfg_no_mfe = ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, use_mfe=False)
        score_no_mfe = calculate_composite_score(
            cai=1.0, gc=60.0, sequence=long_seq, config=cfg_no_mfe
        )
        assert score == pytest.approx(score_no_mfe)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
