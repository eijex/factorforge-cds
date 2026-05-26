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
