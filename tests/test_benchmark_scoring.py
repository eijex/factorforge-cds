import dataclasses
import pandas as pd
import pytest
from benchmarks.config import load_benchmark_config
from benchmarks.scoring import score_cds, canonical_multi_constraint_pass

CFG = load_benchmark_config()


def test_score_cds_full_schema():
    row = score_cds("random_synonymous", "baseline", "p1", "MKT", "ATGAAAACC", CFG, runtime_seconds=0.01)
    expected = {"method","method_type","sequence_id","aa_identity","internal_stop_count",
                "invalid_codon_count","length_multiple_of_three","cai","gc_percent",
                "gc_in_target_range","forbidden_type_iis_site_count","biological_pass",
                "assembly_pass","multi_constraint_pass","runtime_seconds"}
    assert expected.issubset(row.keys())
    assert row["biological_pass"] is True
    assert row["aa_identity"] == 1.0


# --- scoring_contract v1.1 semantic invariant tests ---

def test_multi_constraint_pass_requires_gc_in_range():
    """multi_constraint_pass=True implies gc_in_target_range=True (scoring_contract v1.1)."""
    # Use a very narrow GC window that forces gc_in_target_range=False for a typical CDS
    cfg_narrow = dataclasses.replace(CFG, gc_min=99.0, gc_max=100.0)
    row = score_cds("greedy_cai", "baseline", "t1", "MKT", "ATGAAAACC", cfg_narrow, runtime_seconds=0.01)
    assert row["gc_in_target_range"] is False
    assert row["multi_constraint_pass"] is False, (
        "multi_constraint_pass must be False when gc_in_target_range is False (scoring_contract v1.1)"
    )


def test_multi_constraint_pass_true_only_when_all_primitives_pass():
    """multi_constraint_pass=True only when biological_pass AND assembly_pass AND gc_in_target_range."""
    # MKT codon ATGAAAACC: GC = (0+0+0+0+0+1+0+1+1)/9 = 3/9 ~ 33%
    # Default gc_min=55, gc_max=65 → gc_in_target_range=False → multi_constraint_pass=False
    row = score_cds("greedy_cai", "baseline", "t2", "MKT", "ATGAAAACC", CFG, runtime_seconds=0.01)
    if not row["gc_in_target_range"]:
        assert row["multi_constraint_pass"] is False
    if not row["assembly_pass"]:
        assert row["multi_constraint_pass"] is False


def test_canonical_helper_gc_fail_implies_false():
    """canonical_multi_constraint_pass: gc_in_target_range=False → result=False."""
    df = pd.DataFrame([{
        "biological_pass": True,
        "assembly_pass": True,
        "gc_in_target_range": False,
    }])
    result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)
    assert bool(result.iloc[0]) is False


def test_canonical_helper_all_pass_implies_true():
    """canonical_multi_constraint_pass: all primitives True → result=True."""
    df = pd.DataFrame([{
        "biological_pass": True,
        "assembly_pass": True,
        "gc_in_target_range": True,
    }])
    result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)
    assert bool(result.iloc[0]) is True


def test_canonical_helper_assembly_fail_implies_false():
    """canonical_multi_constraint_pass: assembly_pass=False → result=False."""
    df = pd.DataFrame([{
        "biological_pass": True,
        "assembly_pass": False,
        "gc_in_target_range": True,
    }])
    result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)
    assert bool(result.iloc[0]) is False


def test_canonical_helper_biological_fail_implies_false():
    """canonical_multi_constraint_pass: biological_pass=False → result=False."""
    df = pd.DataFrame([{
        "biological_pass": False,
        "assembly_pass": True,
        "gc_in_target_range": True,
    }])
    result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)
    assert bool(result.iloc[0]) is False


def test_canonical_helper_nan_treated_as_false():
    """canonical_multi_constraint_pass: NaN in any primitive → result=False."""
    import numpy as np
    df = pd.DataFrame([{
        "biological_pass": True,
        "assembly_pass": True,
        "gc_in_target_range": float("nan"),
    }])
    result = canonical_multi_constraint_pass(df, gc_min=55.0, gc_max=65.0)
    assert bool(result.iloc[0]) is False
