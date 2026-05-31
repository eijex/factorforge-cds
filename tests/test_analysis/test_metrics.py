"""Tests for shared sequence-analysis metrics."""

from __future__ import annotations

import pytest

from factorforge.analysis.metrics import (
    HOMOPOLYMER_EXPRESSION_WARN_NT,
    HOMOPOLYMER_SYNTHESIS_WARN_NT,
    amino_acid_identity,
    calculate_cai,
    calculate_first_region_gc,
    calculate_gc,
    calculate_gc_windows,
    codon_usage_profile,
    count_internal_stops,
    detect_forbidden_motifs,
    detect_homopolymers,
    detect_invalid_codons,
    detect_repeats,
    translate_dna,
)


def test_global_gc() -> None:
    assert calculate_gc("") == 0.0
    assert calculate_gc("ATGC") == 50.0
    assert calculate_gc("GGCC") == 100.0


def test_gc_windows_and_first_region_gc() -> None:
    windows = calculate_gc_windows("A" * 30 + "G" * 30, window_size=30, step=30)
    assert [window["gc"] for window in windows] == [0.0, 100.0]

    first = calculate_first_region_gc("G" * 30 + "A" * 60)
    assert first["first_30nt_gc"] == 100.0
    assert first["first_60nt_gc"] == 50.0
    assert first["first_90nt_gc"] == pytest.approx(33.3333333333)


def test_translation_and_amino_acid_identity() -> None:
    assert translate_dna("ATGGCTTTC") == "MAF"
    assert amino_acid_identity("MAF", "ATGGCTTTC") == 1.0
    assert amino_acid_identity("MAF", "ATGGCTTAA") == pytest.approx(2 / 3)


def test_cai_known_weights() -> None:
    weights = {"ATG": 1.0, "GCC": 1.0, "GCT": 0.25}
    assert calculate_cai("ATGGCC", weights) == pytest.approx(1.0)
    assert calculate_cai("ATGGCT", weights) == pytest.approx(0.5)
    assert calculate_cai("ATGG", weights) == 0.0


def test_invalid_codons_and_internal_stops() -> None:
    assert count_internal_stops("ATGTAAGCC") == 1
    invalid = detect_invalid_codons("ATGNN")
    assert invalid[0]["reason"] == "partial"
    assert detect_invalid_codons("ATGNNN")[0]["reason"] == "invalid_base"


def test_homopolymers_repeats_and_forbidden_motifs() -> None:
    assert detect_homopolymers("ATGAAAAAAT", max_run=6)[0]["base"] == "A"
    assert detect_repeats("ATATATGCGC")
    assert detect_forbidden_motifs("ATGGGTCTCTAA", ["GGTCTC"])[0]["motif"] == "GGTCTC"


def test_homopolymer_context_expression_stability() -> None:
    """6 nt run triggers expression warning; output carries context metadata."""
    run_6 = "GCC" + "A" * 6 + "GCC"
    findings = detect_homopolymers(run_6)
    assert len(findings) == 1
    assert findings[0]["context"] == "expression_stability"
    assert findings[0]["threshold_nt"] == HOMOPOLYMER_EXPRESSION_WARN_NT
    assert findings[0]["length"] == 6


def test_homopolymer_thresholds_are_distinct() -> None:
    """Constants confirm the two thresholds serve different purposes."""
    assert HOMOPOLYMER_EXPRESSION_WARN_NT == 6
    assert HOMOPOLYMER_SYNTHESIS_WARN_NT == 8
    assert HOMOPOLYMER_EXPRESSION_WARN_NT < HOMOPOLYMER_SYNTHESIS_WARN_NT


def test_homopolymer_expression_does_not_flag_below_threshold() -> None:
    """5 nt run is below expression threshold — no findings."""
    run_5 = "GCC" + "A" * 5 + "GCC"
    assert detect_homopolymers(run_5) == []


def test_codon_usage_profile() -> None:
    profile = codon_usage_profile("ATGGCCGCC")
    assert profile["GCC"]["count"] == 2
    assert profile["GCC"]["frequency"] == pytest.approx(2 / 3)
    assert profile["ATG"]["aa"] == "M"
