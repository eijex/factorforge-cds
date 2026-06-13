"""Tests for experimental insect host support: Sf9 (S. frugiperda) and Tni (T. ni), issue #23.

Codon tables are PLACEHOLDER data; tests validate structure and wiring, not biological accuracy.
"""

from __future__ import annotations

import pytest

from factorforge.analysis.metrics import translate_dna
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.profile.rules.reverse_translator import ReverseTranslator
from factorforge.engines.profile.utils import get_data_path, load_codon_table
from factorforge.registry.registry_loader import load_registry, resolve_ref


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"

NEW_INSECT_HOSTS = [
    ("sf9", "spodoptera_frugiperda", "Spodoptera frugiperda"),
    ("tni", "trichoplusia_ni", "Trichoplusia ni"),
]


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_INSECT_HOSTS)
def test_codon_table_file_exists(cli_alias, file_stem, expected_organism):
    data_path = get_data_path()
    json_path = data_path / f"{file_stem}_codons.json"
    assert json_path.exists(), f"Missing codon table: {json_path}"


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_INSECT_HOSTS)
def test_codon_table_loads_and_has_64_codons(cli_alias, file_stem, expected_organism):
    payload = load_codon_table(file_stem, get_data_path())
    assert payload["organism"] == expected_organism
    assert len(payload["codons"]) == 64


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_INSECT_HOSTS)
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


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_INSECT_HOSTS)
def test_optimizer_produces_valid_sequence_for_insect_hosts(cli_alias, file_stem, expected_organism):
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
        ReverseTranslator(host="not_an_insect_host")


def test_all_new_insect_hosts_in_registry():
    registry = load_registry()
    hosts = resolve_ref(registry, "parameters.host_profiles")
    for cli_alias, _file_stem, _organism in NEW_INSECT_HOSTS:
        assert cli_alias in hosts, f"Host {cli_alias!r} missing from registry host_profiles"


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_INSECT_HOSTS)
def test_new_insect_hosts_are_experimental_in_registry(cli_alias, file_stem, expected_organism):
    registry = load_registry()
    hosts = resolve_ref(registry, "parameters.host_profiles")
    host = hosts[cli_alias]
    assert host["status"] == "experimental", (
        f"{cli_alias}: expected status='experimental', got {host['status']!r}"
    )
    assert host["scientific_name"] == expected_organism


@pytest.mark.parametrize("cli_alias,file_stem,expected_organism", NEW_INSECT_HOSTS)
def test_insect_codon_tables_are_placeholder(cli_alias, file_stem, expected_organism):
    payload = load_codon_table(file_stem, get_data_path())
    assert payload["total_cds"] == 0, (
        f"{file_stem}: expected total_cds=0 for placeholder table, got {payload['total_cds']}"
    )
    assert payload["total_codons"] == 0, (
        f"{file_stem}: expected total_codons=0 for placeholder table, got {payload['total_codons']}"
    )
    assert "PLACEHOLDER" in payload["source"], (
        f"{file_stem}: expected 'PLACEHOLDER' in source field for unvalidated data"
    )


def test_sf9_preferred_arg_codon_is_aga():
    payload = load_codon_table("spodoptera_frugiperda", get_data_path())
    assert payload["amino_acids"]["R"]["preferred"] == "AGA", (
        "Sf9: AGA is the hallmark insect Arg codon and should be preferred"
    )


def test_tni_preferred_arg_codon_is_aga():
    payload = load_codon_table("trichoplusia_ni", get_data_path())
    assert payload["amino_acids"]["R"]["preferred"] == "AGA", (
        "Tni: AGA is the hallmark insect Arg codon and should be preferred"
    )


def test_tni_is_more_at_biased_than_sf9():
    sf9 = load_codon_table("spodoptera_frugiperda", get_data_path())
    tni = load_codon_table("trichoplusia_ni", get_data_path())
    assert tni["gc_content"]["overall"] < sf9["gc_content"]["overall"], (
        "T. ni should have lower GC content than S. frugiperda (stronger AT-bias)"
    )


def test_sf9_and_tni_both_have_amino_acids_map():
    for file_stem in ("spodoptera_frugiperda", "trichoplusia_ni"):
        payload = load_codon_table(file_stem, get_data_path())
        assert "amino_acids" in payload
        assert len(payload["amino_acids"]) == 21  # 20 aa + stop
