# factorforge/tests/test_quantized_scoring.py
"""Tests for Canonical Fixed-Point Quantization and Deterministic Tie-Breaking."""

import pytest
from factorforge.scoring.quantization import CanonicalQuantizer, PathRank


def test_canonical_quantizer_basic_scale() -> None:
    quantizer = CanonicalQuantizer(scale=1_000_000)
    
    # 1.0 -> 1_000_000
    assert quantizer.quantize(1.0) == 1_000_000
    # 0.834729 -> 834729
    assert quantizer.quantize(0.834729) == 834729
    # Negative log value
    assert quantizer.quantize(-0.123456) == -123456


def test_canonical_quantizer_round_half_even() -> None:
    quantizer = CanonicalQuantizer(scale=10) # scale 10 to test round half even on 0.5 cases
    
    # 0.05 * 10 = 0.5 -> rounds to even integer 0
    assert quantizer.quantize(0.05) == 0
    # 0.15 * 10 = 1.5 -> rounds to even integer 2
    assert quantizer.quantize(0.15) == 2
    # 0.25 * 10 = 2.5 -> rounds to even integer 2
    assert quantizer.quantize(0.25) == 2
    # 0.35 * 10 = 3.5 -> rounds to even integer 4
    assert quantizer.quantize(0.35) == 4


def test_quantize_additive_vector() -> None:
    quantizer = CanonicalQuantizer(scale=1_000_000)
    scores = [0.1, -0.5, 0.999999]
    quantized = quantizer.quantize_additive_vector(scores)
    assert quantized == [100000, -500000, 999999]
    assert sum(quantized) == 599999


def test_path_rank_prefix_comparison() -> None:
    # AAA has lower rank index than AAC
    path_1 = ["ATG", "AAA", "GCT"]
    path_2 = ["ATG", "AAC", "GCT"]
    
    assert PathRank.compare_prefix_paths(path_1, path_2) == -1
    assert PathRank.compare_prefix_paths(path_2, path_1) == 1
    assert PathRank.compare_prefix_paths(path_1, path_1) == 0


def test_path_rank_best_path_selection() -> None:
    # Multiple candidate paths with identical scores must be deterministically chosen by prefix order
    path_a = ["ATG", "CTG", "GCT"]
    path_b = ["ATG", "CTA", "GCT"]  # CTA comes first in alphabetical order (CTA < CTC < CTG < CTT)
    path_c = ["ATG", "CTT", "GCT"]

    candidates = [path_a, path_b, path_c]
    best = PathRank.best_path(candidates)
    assert best == path_b  # CTA comes first lexicographically
