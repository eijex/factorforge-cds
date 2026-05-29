"""Tests for Tobacco BY-2 / N. tabacum host support."""

from __future__ import annotations

from pathlib import Path

import pytest

from factorforge.analysis.metrics import translate_dna
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
