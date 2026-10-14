"""Synthetic regressions for independent checks and privacy-safe exports."""
import pytest
from factorforge.validation.tradeoff_suite import validate_hard_constraints, analyze_edits, compute_tradeoff_comparison

def test_reverse_complement_and_default_motifs():
    assert not validate_hard_constraints("GAGACC", "ET", 0, 1, ["GGTCTC"])
    assert not validate_hard_constraints("GAGACC", "ET", 0, 1)
    assert validate_hard_constraints("GAGACC", "ET", 0, 1, [])

def test_full_construct_junction_and_case():
    assert not validate_hard_constraints("GAG", "E", 0, 1, ["ggtctc"], right_flank="ACC")

def test_unrounded_gc_boundary():
    assert not validate_hard_constraints("GAG", "E", 0, 0.66665, [])
    assert validate_hard_constraints("GAG", "E", 0, 2/3, [])

def test_frame_and_identity():
    assert not validate_hard_constraints("GAGA", "E", 0, 1, [])
    assert not validate_hard_constraints("GAG", "D", 0, 1, [])
    with pytest.raises(ValueError):
        analyze_edits("GAG", "GAGGAG", [1])

def test_export_excludes_sequences():
    result = compute_tradeoff_comparison("MDEF", [2], {"ATG":1,"GAT":1,"GAC":.4,"GAA":1,"GAG":.5,"TTT":1,"TTC":.3}, 0, 1, [])
    assert "cds_1dp" not in result.to_dict()
    assert result.to_dict(include_sequences=True)["cds_1dp"] == result.cds_1dp

def test_invalid_baseline_fails_comparison(monkeypatch):
    from factorforge.validation import tradeoff_suite as suite
    monkeypatch.setattr(suite.DPV2Optimizer, "optimize", lambda *args, **kwargs: {"sequence":"GAG", "score":0})
    result = suite.compute_tradeoff_comparison("D", [1], {"GAT":1,"GAC":.4}, 0, 1, [])
    assert not result.hard_constraints_pass
    assert "FAIL" in result.dashboard_table
