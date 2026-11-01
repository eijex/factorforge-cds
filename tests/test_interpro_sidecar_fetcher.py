"""Unit tests for InterPro Automated Sidecar Fetcher & Domain Resolver."""

from unittest.mock import patch

from factorforge.structure.interpro_sidecar_fetcher import (
    InterProSidecarFetcher,
)
from factorforge.nbfold.annotations import _parse_interpro_sidecar

SAMPLE_MOCK_INTERPRO_RECORD = {
    "results": [
        {
            "metadata": {
                "accession": "IPR000452",
                "name": "Tetraspanin, conserved site",
                "source_database": "interpro",
                "type": "domain",
            },
            "proteins": [
                {
                    "accession": "p21926",
                    "protein_length": 226,
                    "entry_protein_locations": [
                        {
                            "fragments": [
                                {
                                    "start": 113,
                                    "end": 195,
                                }
                            ]
                        }
                    ],
                }
            ],
        },
        {
            "metadata": {
                "accession": "PF00335",
                "name": "Tetraspanin family",
                "source_database": "pfam",
                "type": "family",
            },
            "proteins": [
                {
                    "accession": "p21926",
                    "protein_length": 226,
                    "entry_protein_locations": [
                        {
                            "fragments": [
                                {
                                    "start": 12,
                                    "end": 208,
                                }
                            ]
                        }
                    ],
                }
            ],
        },
    ]
}

SAMPLE_PROTEIN_226 = "M" * 226


def test_interpro_fetcher_with_mock() -> None:
    fetcher = InterProSidecarFetcher()
    with patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value=SAMPLE_MOCK_INTERPRO_RECORD,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P21926", target_protein_sequence=SAMPLE_PROTEIN_226)

    assert sidecar["status"] == "AVAILABLE"
    assert sidecar["reviewed"] is False  # Mandatory human review invariant
    assert sidecar["accession"] == "P21926"
    assert sidecar["protein_length"] == 226
    assert sidecar["biological_validation"] == "NOT_AVAILABLE"

    # Both entries (domain and family) are preserved in entries list for annotation and contract compliance
    assert len(sidecar["entries"]) == 2
    types = [e["entry_type"] for e in sidecar["entries"]]
    assert "domain" in types
    assert "family" in types

    # Invariant: ONLY 'domain' entries generate candidate_regions (interpro_domain)
    # 'family' entries must NOT be treated as domain folding pause boundaries
    assert len(sidecar["candidate_regions"]) == 1
    assert sidecar["candidate_regions"][0]["region_type"] == "interpro_domain"
    assert sidecar["candidate_regions"][0]["start_1based"] == 113
    assert sidecar["candidate_regions"][0]["end_1based"] == 195

    # Recommended positions must include domain boundary 195, but NOT family span 208
    assert 195 in sidecar["recommended_positions"]
    assert 208 not in sidecar["recommended_positions"]


def test_interpro_sidecar_nbfold_contract_compliance() -> None:
    """Verify that generated InterPro sidecar strictly complies with nbfold annotations contract."""
    fetcher = InterProSidecarFetcher()
    with patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value=SAMPLE_MOCK_INTERPRO_RECORD,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P21926", target_protein_sequence=SAMPLE_PROTEIN_226)

    # Must pass nbfold.annotations._parse_interpro_sidecar without exception
    parsed = _parse_interpro_sidecar(
        interpro_data=sidecar,
        protein_sequence=SAMPLE_PROTEIN_226,
        protein_sha256=sidecar["protein_sha256"],
        interpro_sidecar_sha256=None,
        stage="test_verification",
    )

    assert len(parsed) == 2
    accs = [p["accession"] for p in parsed]
    assert "IPR000452" in accs
    assert "PF00335" in accs


def test_interpro_fetcher_mismatched_target_accession_rejected() -> None:
    """When the mock or API response contains proteins for a different accession, reject it."""
    fetcher = InterProSidecarFetcher()
    with patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value=SAMPLE_MOCK_INTERPRO_RECORD,  # contains p21926
    ):
        # Requesting P99999 must detect accession mismatch and return ERROR_ACCESSION_MISMATCH
        res = fetcher.fetch_and_format_sidecar("P99999", target_protein_sequence=SAMPLE_PROTEIN_226)

    assert res["status"] == "ERROR_ACCESSION_MISMATCH"
    assert "none match requested accession" in res["reason"]
    assert len(res["entries"]) == 0
    assert len(res["candidate_regions"]) == 0


def test_interpro_fetcher_no_entries_found() -> None:
    """When InterPro response has empty results, status must be NO_ENTRIES_FOUND, not AVAILABLE."""
    fetcher = InterProSidecarFetcher()
    with patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value={"results": []},
    ):
        res = fetcher.fetch_and_format_sidecar("P21926", target_protein_sequence=SAMPLE_PROTEIN_226)

    assert res["status"] == "NO_ENTRIES_FOUND"
    assert "No InterPro domain or family entries found" in res["reason"]


def test_interpro_fetcher_guards() -> None:
    fetcher = InterProSidecarFetcher()

    # Invalid accession
    res1 = fetcher.fetch_and_format_sidecar("", target_protein_sequence=SAMPLE_PROTEIN_226)
    assert res1["status"] == "ERROR_INVALID_ACCESSION"

    # Empty sequence
    res2 = fetcher.fetch_and_format_sidecar("P21926", target_protein_sequence="   ")
    assert res2["status"] == "ERROR_INVALID_SEQUENCE"


def test_interpro_family_and_repeat_not_treated_as_domain() -> None:
    """Explicitly verify that family, repeat, and active_site entries do not become domain pause candidates."""
    mock_mixed = {
        "results": [
            {
                "metadata": {"accession": "PF00001", "name": "FamA", "type": "family"},
                "proteins": [{"accession": "p21926", "entry_protein_locations": [{"fragments": [{"start": 1, "end": 200}]}]}],
            },
            {
                "metadata": {"accession": "PF00002", "name": "RepB", "type": "repeat"},
                "proteins": [{"accession": "p21926", "entry_protein_locations": [{"fragments": [{"start": 20, "end": 40}]}]}],
            },
            {
                "metadata": {"accession": "PF00003", "name": "DomC", "type": "domain"},
                "proteins": [{"accession": "p21926", "entry_protein_locations": [{"fragments": [{"start": 50, "end": 90}]}]}],
            },
        ]
    }

    fetcher = InterProSidecarFetcher()
    with patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value=mock_mixed,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P21926", target_protein_sequence=SAMPLE_PROTEIN_226)

    assert sidecar["status"] == "AVAILABLE"
    assert len(sidecar["entries"]) == 3
    # Only DomC (domain) must be in candidate_regions
    assert len(sidecar["candidate_regions"]) == 1
    assert sidecar["candidate_regions"][0]["region_type"] == "interpro_domain"
    assert sidecar["candidate_regions"][0]["start_1based"] == 50
    assert sidecar["candidate_regions"][0]["end_1based"] == 90
    assert sidecar["recommended_positions"] == [90]
