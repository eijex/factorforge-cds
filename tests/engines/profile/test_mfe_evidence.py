"""MFE provenance/evidence tests (job 075, 016 audit blocker).

When MFE is not computed (e.g. ViennaRNA unavailable, as on Vercel), downstream
artifacts must report it honestly — never as a misleading 0.0.
"""

from __future__ import annotations

import factorforge.engines.profile.scoring as scoring
from factorforge.engines.profile.scoring import (
    ScoringConfig,
    calculate_composite_score,
    compute_mfe_evidence,
)
from factorforge.engines.profile.pipeline import PipelineResult

SEQ = "ATG" + "GCT" * 20 + "TAA"


def test_evidence_not_computed_when_vienna_unavailable(monkeypatch) -> None:
    monkeypatch.setattr(scoring, "_check_vienna_available", lambda: False)
    ev = compute_mfe_evidence(SEQ, profile="balanced")

    assert ev["mfe_kcal_mol"] is None
    assert ev["mfe_status"] == "not_computed"
    assert ev["mfe_used"] is False
    assert "ViennaRNA" in ev["mfe_warning"]
    assert ev["score_components"]["mfe_used"] is False
    assert ev["score_components"]["cai_used"] is True


def test_evidence_disabled_profile(monkeypatch) -> None:
    # Even if ViennaRNA is available, a use_mfe=False config must report honestly.
    monkeypatch.setattr(scoring, "_check_vienna_available", lambda: True)
    cfg = ScoringConfig(w_cai=0.6, w_gc=0.4, use_mfe=False)
    ev = compute_mfe_evidence(SEQ, config=cfg)

    assert ev["mfe_status"] == "not_computed"
    assert ev["mfe_used"] is False
    assert "disabled" in ev["mfe_warning"]


def test_evidence_no_sequence() -> None:
    ev = compute_mfe_evidence(None, profile="balanced")
    assert ev["mfe_kcal_mol"] is None
    assert ev["mfe_used"] is False


def test_composite_score_unchanged_without_mfe(monkeypatch) -> None:
    # Score value must not depend on the evidence helper; MFE-absent path is
    # handled inside calculate_composite_score and stays a plain float.
    monkeypatch.setattr(scoring, "_check_vienna_available", lambda: False)
    score = calculate_composite_score(cai=0.85, gc=60.0, sequence=SEQ, profile="balanced")
    assert isinstance(score, float)
    assert 0.0 <= score <= 1.0


def test_export_features_reports_null_not_zero_when_not_computed() -> None:
    # Design Package must show null (not 0.0) for an uncomputed MFE.
    result = PipelineResult(
        sequence=SEQ,
        metadata={
            "construct_id": "CF-TEST",
            "profile": "balanced",
            "metrics": {
                "cai": 0.85,
                "gc": 60.0,
                "mfe_kcal_mol": None,
                "mfe_status": "not_computed",
                "mfe_used": False,
            },
        },
    )
    features = result.export_features()
    assert features["mfe_kcal_mol"] is None
    assert features["mfe_status"] == "not_computed"
    assert features["mfe_used"] is False


def test_export_features_reports_value_when_computed() -> None:
    result = PipelineResult(
        sequence=SEQ,
        metadata={
            "construct_id": "CF-TEST",
            "profile": "balanced",
            "metrics": {
                "cai": 0.85,
                "gc": 60.0,
                "mfe_kcal_mol": -12.345,
                "mfe_status": "computed",
                "mfe_used": True,
            },
        },
    )
    features = result.export_features()
    assert features["mfe_kcal_mol"] == -12.35
    assert features["mfe_status"] == "computed"
    assert features["mfe_used"] is True
