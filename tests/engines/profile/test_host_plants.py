"""Tests for new plant host support: Arabidopsis, Tomato, Lemna, Wolffia (issue #24)."""

from __future__ import annotations

import pytest

from factorforge.analysis.metrics import translate_dna
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.engines.profile.utils import get_data_path, load_codon_table
from factorforge.registry.registry_loader import load_registry, resolve_ref


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"

NEW_PLANT_HOSTS = [
    ("arabidopsis", "arabidopsis_thaliana", "Arabidopsis thaliana"),
    ("tomato", "solanum_lycopersicum", "Solanum lycopersicum"),
    ("lemna", "lemna_minor", "Lemna minor"),
    ("wolffia", "wolffia_globosa", "Wolffia globosa"),
]


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_PLANT_HOSTS)
def test_codon_table_file_exists(cli_alias, file_stem, expected_organism):
    data_path = get_data_path()
    json_path = data_path / f"{file_stem}_codons.json"
    assert json_path.exists(), f"Missing codon table: {json_path}"


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_PLANT_HOSTS)
def test_codon_table_loads_and_has_64_codons(cli_alias, file_stem, expected_organism):
    payload = load_codon_table(file_stem, get_data_path())
    assert payload["organism"] == expected_organism
    assert len(payload["codons"]) == 64


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_PLANT_HOSTS)
def test_codon_frequencies_normalized(cli_alias, file_stem, expected_organism):
    payload = load_codon_table(file_stem, get_data_path())
    codons = payload["codons"]

    aa_to_freqs: dict[str, list[float]] = {}
    for codon_data in codons.values():
        aa = codon_data["aa"]
        aa_to_freqs.setdefault(aa, []).append(codon_data["frequency"])

    for aa, freqs in aa_to_freqs.items():
        total = sum(freqs)
        assert abs(total - 1.0) < 0.01, (
            f"{file_stem}: frequencies for aa={aa!r} sum to {total:.4f}, expected ~1.0"
        )


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_PLANT_HOSTS)
def test_optimizer_produces_valid_sequence_for_new_hosts(cli_alias, file_stem, expected_organism):
    optimizer = RuleBasedOptimizer()
    result = optimizer.optimize(
        SAMPLE_PROTEIN,
        profile="balanced",
        host=file_stem,
        scan_mode="fast",
    )

    assert result.sequence.startswith("ATG")
    assert len(result.sequence) == len(SAMPLE_PROTEIN) * 3
    assert translate_dna(result.sequence).rstrip("*") == SAMPLE_PROTEIN
    assert result.metadata["host"] == file_stem
    assert "cai" in result.metrics
    assert "gc_percent" in result.metrics


def test_invalid_host_still_raises_file_not_found():
    with pytest.raises(FileNotFoundError):
        ReverseTranslator(host="not_a_real_host")


def test_all_new_plant_hosts_in_registry():
    registry = load_registry()
    hosts = resolve_ref(registry, "parameters.host_profiles")
    for cli_alias, _file_stem, _organism in NEW_PLANT_HOSTS:
        assert cli_alias in hosts, f"Host {cli_alias!r} missing from registry host_profiles"


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_PLANT_HOSTS)
def test_new_plant_hosts_are_experimental_in_registry(cli_alias, file_stem, expected_organism):
    registry = load_registry()
    hosts = resolve_ref(registry, "parameters.host_profiles")
    host = hosts[cli_alias]
    assert host["status"] == "experimental", (
        f"{cli_alias}: expected status='experimental', got {host['status']!r}"
    )
    assert host["scientific_name"] == expected_organism


def test_arabidopsis_preferred_codon_is_a_t_biased():
    payload = load_codon_table("arabidopsis_thaliana", get_data_path())
    aa = payload["amino_acids"]
    assert aa["L"]["preferred"] == "CTT"
    assert aa["A"]["preferred"] == "GCT"
    assert aa["P"]["preferred"] == "CCT"


def test_tomato_preferred_codons_are_gc_biased():
    payload = load_codon_table("solanum_lycopersicum", get_data_path())
    aa = payload["amino_acids"]
    assert aa["L"]["preferred"] == "CTG"
    assert aa["A"]["preferred"] == "GCC"
    assert aa["S"]["preferred"] == "AGC"


def test_lemna_preferred_codon_leu_is_ctc():
    payload = load_codon_table("lemna_minor", get_data_path())
    assert payload["amino_acids"]["L"]["preferred"] == "CTC"


def test_wolffia_codon_table_has_expected_organism_field():
    payload = load_codon_table("wolffia_globosa", get_data_path())
    assert payload["organism"] == "Wolffia globosa"
    assert len(payload["codons"]) == 64
