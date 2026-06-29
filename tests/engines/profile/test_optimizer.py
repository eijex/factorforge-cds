"""
Unit tests for RuleBasedOptimizer
"""

import sys
from pathlib import Path

import pytest

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.optimizer import RuleBasedOptimizer


@pytest.fixture
def optimizer():
    """Create optimizer instance"""
    return RuleBasedOptimizer()


@pytest.fixture
def sample_protein():
    """Short protein sequence for tests"""
    return "MAKL"


@pytest.mark.parametrize(
    "profile",
    ["balanced", "high_cai", "gc_target", "assembly_friendly"],
)
def test_optimize_profiles(optimizer, sample_protein, profile):
    """Optimize with all supported profiles"""
    result = optimizer.optimize(sample_protein, profile=profile)

    assert result.sequence
    assert len(result.sequence) == len(sample_protein) * 3
    assert "cai" in result.metrics
    assert "gc_content" in result.metrics
    assert "score" in result.metrics


@pytest.mark.parametrize(
    ("protein", "target_gc_min", "target_gc_max", "seed"),
    [
        ("MAKL", 55.0, 65.0, 174),
        ("MSKGEELFTGVVPILVELD", 40.0, 47.0, 174),
        ("GGGGKKKK", 30.0, 70.0, 174),
    ],
)
def test_balanced_gc_target_reached_matches_achieved_gc(
    optimizer, protein, target_gc_min, target_gc_max, seed
):
    """Balanced GC target flag is a direct observation of achieved GC%.

    This is a fixed-fixture logic test, not a sampling-rate estimate.
    """
    result = optimizer.optimize(
        protein,
        profile="balanced",
        target_gc_min=target_gc_min,
        target_gc_max=target_gc_max,
        seed=seed,
        scan_mode="fast",
    )
    achieved_gc_percent = result.metrics["gc_percent"]

    assert result.metrics["gc_content"] == achieved_gc_percent
    assert result.metrics["requested_gc_min_percent"] == target_gc_min
    assert result.metrics["requested_gc_max_percent"] == target_gc_max
    assert result.metrics["gc_target_reached"] == (
        target_gc_min <= achieved_gc_percent <= target_gc_max
    )


@pytest.mark.parametrize("profile", ["high_cai", "gc_target", "assembly_friendly"])
def test_gc_target_reached_metrics_are_balanced_only(optimizer, sample_protein, profile):
    """GC target attainment fields are scoped to balanced profile output."""
    result = optimizer.optimize(sample_protein, profile=profile, seed=174, scan_mode="fast")

    assert "gc_target_reached" not in result.metrics
    assert "requested_gc_min_percent" not in result.metrics
    assert "requested_gc_max_percent" not in result.metrics


def test_optimize_invalid_input(optimizer):
    """Reject invalid input"""
    with pytest.raises(ValueError):
        optimizer.optimize("12345", profile="balanced")


def test_optimize_invalid_profile_raises(optimizer, sample_protein):
    """Unknown profile should raise instead of silently falling back."""
    with pytest.raises(ValueError, match="Unknown profile"):
        optimizer.optimize(sample_protein, profile="fast")


def test_optimize_fast_scan_mode(optimizer, sample_protein):
    """Fast scan mode should skip heavier scanners."""
    result = optimizer.optimize(sample_protein, profile="balanced", scan_mode="fast")
    scan_results = result.metadata["scan_results"]

    assert result.metadata["scan_mode"] == "fast"
    assert "polya" in scan_results
    assert "repeats" not in scan_results
    assert "dinucleotides" not in scan_results


def test_optimize_batch(optimizer, sample_protein):
    """Batch optimization should preserve order and include input ids."""
    batch = [
        {"id": "a", "sequence": sample_protein},
        {"id": "b", "sequence": sample_protein},
    ]
    results = optimizer.optimize_batch(batch, profile="balanced", scan_mode="fast")

    assert len(results) == 2
    assert results[0].metadata["input_id"] == "a"
    assert results[1].metadata["input_id"] == "b"
    assert len(results[0].sequence) == len(sample_protein) * 3


@pytest.mark.parametrize("seq,expected_len", [
    ("MSCSNYRRC", 27),   # 9 AA — chars overlap with IUPAC ambiguous DNA
    ("MMCWCKMMMS", 30),  # 10 AA
    ("MDSARNNKD", 27),   # 9 AA
])
def test_optimize_short_iupac_overlap_proteins(optimizer, seq, expected_len):
    """Proteins whose AA codes overlap with IUPAC ambiguous DNA must produce a valid CDS.

    Before the fix these sequences were misclassified as DNA and returned unchanged,
    producing invalid_codon_count > 0 and aa_identity = 0.0 in the benchmark.
    (see test_input_validator.py for the benchmark case this regression-guards)
    """
    result = optimizer.optimize(seq, profile="balanced")
    assert len(result.sequence) == expected_len, (
        f"Expected CDS length {expected_len} for {seq}, got {len(result.sequence)}"
    )
