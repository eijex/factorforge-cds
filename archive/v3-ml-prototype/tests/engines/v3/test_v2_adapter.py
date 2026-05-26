"""Tests for the v3-alpha v2 adapter."""

from __future__ import annotations

import pytest

from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2
from fixtures.benchmark_proteins import BENCHMARK_PROTEINS


def test_v2_adapter_returns_dna_sequence() -> None:
    result = optimize_with_v2("MAST")

    assert result["engine"] == "v2"
    assert result["dna_sequence"]
    assert len(result["dna_sequence"]) == len("MAST") * 3


def test_v2_adapter_output_length_divisible_by_three() -> None:
    result = optimize_with_v2(BENCHMARK_PROTEINS["reporter_like"])

    assert len(result["dna_sequence"]) % 3 == 0


def test_v2_adapter_preserves_amino_acids_and_metrics() -> None:
    result = optimize_with_v2(BENCHMARK_PROTEINS["short_synthetic"])

    assert result["validator"]["amino_acid_identity"] == 1.0
    assert result["validator"]["passed"] is True
    assert "cai" in result["metrics"]
    assert "gc" in result["metrics"]
    assert result["validator"]["internal_stop_count"] == 0


def test_v2_adapter_reports_failures_cleanly() -> None:
    with pytest.raises(Exception) as exc_info:
        optimize_with_v2("MZX")

    assert "Invalid amino acid" in str(exc_info.value) or "No valid codon candidates" in str(
        exc_info.value
    )
