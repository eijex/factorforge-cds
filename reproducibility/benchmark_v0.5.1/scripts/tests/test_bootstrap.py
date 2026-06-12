import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "analysis"))
from bootstrap_ci import build_sequence_pass_matrix, paired_bootstrap


def _make_df(n_seq=100, n_replicates_random=3, seed=42):
    rng = np.random.default_rng(seed)
    rows = []
    for m in ["factorforge_assembly_friendly", "greedy_cai", "random_synonymous"]:
        reps = n_replicates_random if m == "random_synonymous" else 1
        for seq_i in range(n_seq):
            for rep in range(1, reps + 1):
                pass_val = rng.random() < (0.65 if "assembly" in m else 0.3)
                rows.append({
                    "method": m,
                    "sequence_id": f"seq{seq_i:04d}",
                    "replicate": rep,
                    "multi_constraint_pass": pass_val,
                })
    return pd.DataFrame(rows)


def test_build_sequence_pass_matrix_shapes():
    df = _make_df(n_seq=50)
    seq_pass = build_sequence_pass_matrix(df)
    assert "factorforge_assembly_friendly" in seq_pass
    assert "random_synonymous" in seq_pass
    assert seq_pass["factorforge_assembly_friendly"].shape == (50,)
    assert seq_pass["random_synonymous"].shape == (50,)


def test_random_synonymous_aggregated_between_0_and_1():
    df = _make_df(n_seq=50)
    seq_pass = build_sequence_pass_matrix(df)
    arr = seq_pass["random_synonymous"]
    assert arr.min() >= 0.0
    assert arr.max() <= 1.0
    # Should have non-integer values (averages of 3 booleans)
    has_fractional = any(v not in (0.0, 1.0) for v in arr)
    assert has_fractional, "Expected fractional values from replicate averaging"


def test_paired_bootstrap_metadata():
    df = _make_df(n_seq=50)
    seq_pass = build_sequence_pass_matrix(df)
    result = paired_bootstrap(seq_pass, n_resamples=100, seed=320)
    assert result["bootstrap_unit"] == "sequence_id"
    assert result["paired_across_methods"] is True
    assert result["n_resamples"] == 100
    assert result["random_synonymous_aggregation"] == "mean_across_replicates_per_sequence"


def test_paired_bootstrap_ci_bounds():
    df = _make_df(n_seq=200, seed=99)
    seq_pass = build_sequence_pass_matrix(df)
    result = paired_bootstrap(seq_pass, n_resamples=200, seed=320)
    for m, ci in result["per_method"].items():
        assert ci["ci_lower_2_5"] <= ci["observed_rate"] <= ci["ci_upper_97_5"], \
            f"{m}: observed rate outside CI"


def test_difference_cis_present():
    df = _make_df(n_seq=100)
    seq_pass = build_sequence_pass_matrix(df)
    result = paired_bootstrap(seq_pass, n_resamples=100, seed=320)
    assert "factorforge_assembly_friendly - greedy_cai" in result["difference_cis"]
    assert "factorforge_assembly_friendly - random_synonymous" in result["difference_cis"]
