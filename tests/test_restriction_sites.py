"""Tests for custom restriction-site detection and synonymous removal."""

from __future__ import annotations

import pytest

from api.optimize import handler
from factorforge.engines.v2.rules.domesticator import Domesticator
from factorforge.engines.v2.utils import get_data_path, load_codon_table
from factorforge.ml.metrics import translate_dna
from factorforge.utils.restriction_sites import (
    detect_restriction_sites,
    domesticate_custom_sites,
)


@pytest.fixture(scope="module")
def codon_table() -> dict:
    return load_codon_table("nbenthamiana", get_data_path())


def test_ecori_removal_feasible(codon_table: dict) -> None:
    seq = "ATGGAATTCGAGTAA"
    original_aa = translate_dna(seq).rstrip("*")

    result = domesticate_custom_sites(
        seq,
        [{"name": "EcoRI", "sequence": "GAATTC"}],
        codon_table,
    )

    assert result["removed"] == [{"name": "EcoRI", "position": 3, "substitution": "GAA->GAG"}]
    assert result["unresolved"] == []
    assert (
        detect_restriction_sites(result["sequence"], [{"name": "EcoRI", "sequence": "GAATTC"}])
        == []
    )
    assert translate_dna(result["sequence"]).rstrip("*") == original_aa


def test_xhoi_removal_feasible(codon_table: dict) -> None:
    seq = "ATGCTCGAGTAA"
    original_aa = translate_dna(seq).rstrip("*")

    result = domesticate_custom_sites(
        seq,
        [{"name": "XhoI", "sequence": "CTCGAG"}],
        codon_table,
    )

    assert result["removed"][0]["name"] == "XhoI"
    assert result["unresolved"] == []
    assert (
        detect_restriction_sites(result["sequence"], [{"name": "XhoI", "sequence": "CTCGAG"}]) == []
    )
    assert translate_dna(result["sequence"]).rstrip("*") == original_aa


def test_reverse_complement_detection() -> None:
    seq = "ATGGAGACCGCC"

    hits = detect_restriction_sites(
        seq,
        [{"name": "BsaI_custom", "sequence": "GGTCTC"}],
    )

    assert hits == [{"name": "BsaI_custom", "sequence": "GGTCTC", "position": 3, "strand": "rc"}]


def test_overlapping_sites(codon_table: dict) -> None:
    seq = "ATGGAATTCGAGTAA"
    sites = [
        {"name": "EcoRI", "sequence": "GAATTC"},
        {"name": "Overlap", "sequence": "AATTCG"},
    ]

    result = domesticate_custom_sites(seq, sites, codon_table)

    assert result["unresolved"] == []
    assert len(result["detected"]) == 2
    assert len(result["removed"]) >= 1
    assert detect_restriction_sites(result["sequence"], sites) == []
    assert translate_dna(result["sequence"]).rstrip("*") == translate_dna(seq).rstrip("*")


def test_removal_introduces_new_site() -> None:
    codon_table = {
        "amino_acids": {
            "G": {"codons": ["GGA", "GGT", "GGC"]},
            "L": {"codons": ["CTC"]},
        },
        "codons": {
            "GGA": {"aa": "G", "frequency": 0.1},
            "GGT": {"aa": "G", "frequency": 0.9},
            "GGC": {"aa": "G", "frequency": 0.8},
            "CTC": {"aa": "L", "frequency": 1.0},
        },
    }
    seq = "GGACTC"

    result = domesticate_custom_sites(
        seq,
        [{"name": "Custom", "sequence": "GGACTC"}],
        codon_table,
    )

    assert result["sequence"] == "GGCCTC"
    assert result["removed"] == [{"name": "Custom", "position": 0, "substitution": "GGA->GGC"}]
    assert (
        detect_restriction_sites(
            result["sequence"],
            [
                {"name": "Custom", "sequence": "GGACTC"},
                {"name": "BsaI", "sequence": "GGTCTC"},
            ],
        )
        == []
    )


def test_unresolved_site() -> None:
    codon_table = {
        "amino_acids": {
            "G": {"codons": ["GGT"]},
            "L": {"codons": ["CTC"]},
        },
        "codons": {
            "GGT": {"aa": "G", "frequency": 1.0},
            "CTC": {"aa": "L", "frequency": 1.0},
        },
    }
    seq = "GGTCTC"

    result = domesticate_custom_sites(
        seq,
        [{"name": "BsaI_custom", "sequence": "GGTCTC"}],
        codon_table,
    )

    assert result["sequence"] == seq
    assert result["removed"] == []
    assert result["unresolved"] == [
        {"name": "BsaI_custom", "position": 0, "reason": "no_synonymous_substitution"}
    ]


def test_cds_boundary_site(codon_table: dict) -> None:
    seq = "ATGGCCTAA"

    result = domesticate_custom_sites(
        seq,
        [{"name": "Boundary", "sequence": "CCTAA"}],
        codon_table,
    )

    assert result["sequence"] == seq
    assert result["removed"] == []
    assert result["unresolved"] == [{"name": "Boundary", "position": 4, "reason": "outside_cds"}]


def test_non_cds_warning_only(codon_table: dict) -> None:
    seq = "ATGGGATGA"

    result = domesticate_custom_sites(
        seq,
        [{"name": "StopBoundary", "sequence": "ATGA"}],
        codon_table,
    )

    assert result["sequence"] == seq
    assert result["removed"] == []
    assert result["unresolved"] == [
        {"name": "StopBoundary", "position": 5, "reason": "outside_cds"}
    ]


def test_assembly_friendly_regression() -> None:
    domesticator = Domesticator()
    result = domesticator.domesticate("ATGGGTCTCGAGGAGCTGTTC", "golden_gate")

    assert "GGTCTC" not in result["domesticated_seq"]
    assert result["removed_sites"][0]["enzyme"] == "BsaI"

    h = object.__new__(handler)
    api_result = h.optimize_sequence(
        "MEFE",
        "gc_target",
        False,
        False,
        False,
        objective=None,
        return_candidates=False,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )

    assert "custom_restriction_sites" not in api_result
