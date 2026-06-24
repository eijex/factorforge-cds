"""Tests for Tobacco BY-2 / N. tabacum host support."""

from __future__ import annotations

from pathlib import Path

import pytest

from factorforge.analysis.metrics import translate_dna
from factorforge.engines.profile.exporter import SequenceExporter
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.engines.profile.utils import get_data_path, load_codon_table


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"


def test_ntabacum_codon_table_loads() -> None:
    payload = load_codon_table("ntabacum", get_data_path())

    assert payload["organism"] == "Nicotiana tabacum"
    assert len(payload["codons"]) == 64
    assert payload["amino_acids"]["F"]["preferred"] == "TTT"


def test_by2_host_optimizer_returns_valid_sequence() -> None:
    optimizer = RuleBasedOptimizer()

    result = optimizer.optimize(
        SAMPLE_PROTEIN,
        profile="balanced",
        host="ntabacum",
        scan_mode="fast",
    )

    assert result.sequence.startswith("ATG")
    assert len(result.sequence) == len(SAMPLE_PROTEIN) * 3
    assert translate_dna(result.sequence).rstrip("*") == SAMPLE_PROTEIN
    assert result.metadata["host"] == "ntabacum"
    assert "cai" in result.metrics
    assert "gc_percent" in result.metrics


def test_invalid_host_raises_file_not_found() -> None:
    with pytest.raises(FileNotFoundError):
        ReverseTranslator(host="not_a_host")


def test_nbenthamiana_and_by2_hosts_both_produce_valid_sequences() -> None:
    optimizer = RuleBasedOptimizer()

    nb_result = optimizer.optimize(
        SAMPLE_PROTEIN,
        profile="balanced",
        host="nbenthamiana",
        scan_mode="fast",
    )
    by2_result = optimizer.optimize(
        SAMPLE_PROTEIN,
        profile="balanced",
        host="ntabacum",
        scan_mode="fast",
    )

    assert translate_dna(nb_result.sequence).rstrip("*") == SAMPLE_PROTEIN
    assert translate_dna(by2_result.sequence).rstrip("*") == SAMPLE_PROTEIN
    assert nb_result.metadata["host"] == "nbenthamiana"
    assert by2_result.metadata["host"] == "ntabacum"
    assert Path(get_data_path(), "ntabacum_codons.json").exists()


WITNESS_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"

# Pinned 2026-06-23 by running the pre-guard code directly (seed=42 fixes
# balanced/assembly_friendly non-determinism, confirmed by running each
# combination twice and comparing). These are not derived by hand. The
# host/strategy compatibility guard below must not change any of these
# values — it only adds rejection/disclosure around the existing engine,
# never touches profile output.
NBENTHAMIANA_BASELINE_SEED42 = {
    "balanced": "ATGAAGACCGCCTACATCGCCAAGCAGCGGCAGATCAGCTTCGTCAAGAGCCATTTCAGCAGGCAGCTGGAAGAGAGGCTGGGACTGATCGAGGTGCAG",
    "gc_target": "ATGAAGACCGCCTACATCGCCAAGCAGCGCCAGATCAGCTTCGTGAAGAGCCACTTCAGCCGCCAGCTGGAGGAGAGGCTGGGACTTATCGAGGTGCAG",
    "assembly_friendly": "ATGAAAACCGCCTACATCGCCAAGCAGAGGCAGATCAGCTTCGTGAAAAGCCACTTCAGCCGGCAGCTCGAGGAGAGGCTGGGCCTCATCGAGGTGCAG",
    "high_cai": "ATGAAAACTGCTTATATTGCTAAACAAAGACAAATTTCTTTCGTTAAATCTCATTTCTCTAGACAACTTGAAGAAAGACTTGGTCTTATTGAAGTTCAA",
}
NTABACUM_BASELINE_SEED42 = {
    "balanced": "ATGAAGACAGCTTATATTGCCAAACAAAGACAAATCTCTTTTGTTAAAAGCCATTTTAGCCGACAACTTGAGGAAAGGCTTGGACTTATTGAAGTCCAA",
    "gc_target": "ATGAAGACCGCCTACATCGCCAAGCAGCGCCAGATCTCCTTCGTGAAGTCCCACTTCTCCCGCCAGCTCGAGGAGAGGCTCGGACTTATCGAGGTGCAG",
    "assembly_friendly": "ATGAAAACAGCTTATATCGCCAAGCAAAGACAAATTTCATTCGTGAAATCTCATTTTAGCAGACAACTTGAGGAACGTCTTGGATTGATTGAGGTTCAA",
    # Intentionally identical to NBENTHAMIANA_BASELINE_SEED42["high_cai"]:
    # high_cai always optimizes against the nbenthamiana-only golden-set
    # reference and has no BY-2 equivalent.
    "high_cai": "ATGAAAACTGCTTATATTGCTAAACAAAGACAAATTTCTTTCGTTAAATCTCATTTCTCTAGACAACTTGAAGAAAGACTTGGTCTTATTGAAGTTCAA",
}


@pytest.mark.parametrize("profile", ["balanced", "gc_target", "assembly_friendly", "high_cai"])
def test_nbenthamiana_baseline_unchanged(profile: str) -> None:
    optimizer = RuleBasedOptimizer()
    result = optimizer.optimize(WITNESS_PROTEIN, profile=profile, host="nbenthamiana", seed=42)
    assert result.sequence == NBENTHAMIANA_BASELINE_SEED42[profile]


@pytest.mark.parametrize("profile", ["balanced", "gc_target", "assembly_friendly", "high_cai"])
def test_ntabacum_baseline_unchanged(profile: str) -> None:
    optimizer = RuleBasedOptimizer()
    result = optimizer.optimize(WITNESS_PROTEIN, profile=profile, host="ntabacum", seed=42)
    assert result.sequence == NTABACUM_BASELINE_SEED42[profile]


@pytest.mark.parametrize("profile", ["balanced", "gc_target", "assembly_friendly"])
def test_host_aware_profiles_differ_by_host(profile: str) -> None:
    """Regression guard for a gap found during the host/strategy guard
    investigation: test_nbenthamiana_and_by2_hosts_both_produce_valid_sequences() never
    asserted the two hosts produce *different* DNA for host-aware profiles,
    which is why high_cai's host-blindness went undetected."""
    assert NBENTHAMIANA_BASELINE_SEED42[profile] != NTABACUM_BASELINE_SEED42[profile]


def test_high_cai_is_host_invariant_by_design() -> None:
    """high_cai is anchored to the nbenthamiana-only golden-set reference
    and has no BY-2 equivalent — this equality is the documented capability
    boundary, not a bug."""
    assert NBENTHAMIANA_BASELINE_SEED42["high_cai"] == NTABACUM_BASELINE_SEED42["high_cai"]


def test_high_cai_logs_warning_for_non_default_host(caplog: pytest.LogCaptureFixture) -> None:
    """Library-direct calls don't reject high_cai+non-default-host (that
    would conflict with the host-invariant-by-design contract above and the
    benchmark's codon_table_path injection path) but should surface a
    warning so silent host-invariance isn't mistaken for host-awareness."""
    optimizer = RuleBasedOptimizer()
    with caplog.at_level("WARNING"):
        optimizer.optimize(WITNESS_PROTEIN, profile="high_cai", host="ntabacum", seed=42)
    assert any("high_cai" in record.message and "ntabacum" in record.message for record in caplog.records)


def test_high_cai_no_warning_for_default_host(caplog: pytest.LogCaptureFixture) -> None:
    optimizer = RuleBasedOptimizer()
    with caplog.at_level("WARNING"):
        optimizer.optimize(WITNESS_PROTEIN, profile="high_cai", host="nbenthamiana", seed=42)
    assert not caplog.records


def test_no_warning_for_non_high_cai_profile_on_non_default_host(
    caplog: pytest.LogCaptureFixture,
) -> None:
    optimizer = RuleBasedOptimizer()
    with caplog.at_level("WARNING"):
        optimizer.optimize(WITNESS_PROTEIN, profile="balanced", host="ntabacum", seed=42)
    assert not caplog.records


def test_by2_genbank_organism_comes_from_feature_registry() -> None:
    exporter = SequenceExporter()

    try:
        genbank = exporter.export_genbank(
            "ATGAAAAAAAAA",
            {
                "gene_name": "by2_test",
                "protein_seq": "MKKK",
                "host_profile": "by2",
                "cai": 0.8,
                "gc": 16.7,
            },
        )
    except ImportError:
        pytest.skip("Biopython not installed")

    assert "Nicotiana tabacum" in genbank
