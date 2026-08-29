import pytest
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.lm.inference import LMEngineAdapter

def test_cross_engine_evaluation_contract():
    """Verify that all engines produce matching EvaluationResult formats."""
    # A tiny test protein
    protein = "MVLSA"
    
    # 1. Profile engine
    profile_engine = RuleBasedOptimizer()
    prof_res = profile_engine.optimize(protein, profile="balanced", host="nbenthamiana", terminal_stop_policy="append")
    
    # 2. LM engine
    lm_engine = LMEngineAdapter()
    lm_res = lm_engine.optimize(protein, host="nbenthamiana", terminal_stop_policy="append")
    
    # Assert structural consistency
    assert "evaluation_report" in prof_res.metadata
    assert "evaluation_report" in lm_res.metadata
    
    prof_eval = prof_res.metadata["evaluation_report"]
    lm_eval = lm_res.metadata["evaluation_report"]
    
    # Both should have exactly the same schema for sequence_integrity
    assert set(prof_eval["sequence_integrity"].keys()) == {"aa_identity", "frame_valid", "internal_stop_count"}
    assert set(lm_eval["sequence_integrity"].keys()) == {"aa_identity", "frame_valid", "internal_stop_count"}
    
    # Both should have metric schema
    assert set(prof_eval["metrics"].keys()) == {"cai", "gc_percent", "mfe"}
    assert set(lm_eval["metrics"].keys()) == {"cai", "gc_percent", "mfe"}
    
    # Both should pass basic invariants (or fail depending on GC bounds for short proteins)
    assert isinstance(prof_res.metadata["validator_passed"], bool)
    assert isinstance(lm_res.metadata["validator_passed"], bool)
