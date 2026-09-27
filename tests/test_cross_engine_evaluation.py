from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.lm.inference import LMEngineAdapter
from factorforge.evaluation.models import EvaluationResult


def test_cross_engine_evaluation_contract():
    """Verify that all engines produce matching EvaluationResult formats."""
    # A tiny test protein
    protein = "MVLSA"

    # 1. Profile engine
    profile_engine = RuleBasedOptimizer()
    prof_res = profile_engine.optimize(
        protein, profile="balanced", host="nbenthamiana", terminal_stop_policy="append"
    )

    # 2. LM engine
    lm_engine = LMEngineAdapter()
    lm_res = lm_engine.optimize(protein, host="nbenthamiana", terminal_stop_policy="append")

    # Assert structural consistency
    assert "evaluation_report" in prof_res.metadata
    assert "evaluation_report" in lm_res.metadata

    prof_eval = prof_res.metadata["evaluation_report"]
    lm_eval = lm_res.metadata["evaluation_report"]

    # Both should have exactly the same schema for sequence_integrity
    assert set(prof_eval["sequence_integrity"].keys()) == {
        "aa_identity",
        "frame_valid",
        "internal_stop_count",
    }
    assert set(lm_eval["sequence_integrity"].keys()) == {
        "aa_identity",
        "frame_valid",
        "internal_stop_count",
    }

    # Both should have metric schema
    metric_fields = {
        "cai",
        "gc_percent",
        "mfe",
        "mfe_status",
        "mfe_reason",
        "mfe_warning",
        "mfe_5p_window",
        "mfe_5p_status",
        "cai_5p_ramp",
        "cai_body",
        "gc_5p_ramp_percent",
        "gc_body_percent",
        "context_digest",
    }
    for result, report in ((prof_res, prof_eval), (lm_res, lm_eval)):
        assert set(report["metrics"]) == metric_fields
        evaluated = EvaluationResult.model_validate(report)
        assert evaluated.sequence_integrity.aa_identity == 1.0
        assert evaluated.sequence_integrity.frame_valid is True
        assert evaluated.sequence_integrity.internal_stop_count == 0
        assert evaluated.metrics.mfe_status in {"computed", "not_computed"}
        if evaluated.metrics.mfe_status == "not_computed":
            assert evaluated.metrics.mfe is None
        assert result.metadata["validator_passed"] is evaluated.passed

    # Both should pass basic invariants (or fail depending on GC bounds for short proteins)
    assert isinstance(prof_res.metadata["validator_passed"], bool)
    assert isinstance(lm_res.metadata["validator_passed"], bool)
