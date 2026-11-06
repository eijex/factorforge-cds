"""Unit tests for UniProt Automated Sidecar Fetcher & Feature Resolver."""

from unittest.mock import patch
from pathlib import Path

from factorforge.structure.uniprot_sidecar_fetcher import (
    UniProtSidecarFetcher,
)
from factorforge.nbfold.annotations import analyze_uniprot_tm_signal_annotations

SAMPLE_MOCK_UNIPROT_RECORD = {
    "primaryAccession": "P99999",
    "entryAudit": {"entryVersion": 10, "sequenceVersion": 2},
    "sequence": {
        "value": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAEEGPKLVASQAALA",
        "length": 606,
    },
    "features": [
        {
            "type": "Signal",
            "location": {
                "start": {"value": 1, "modifier": "EXACT"},
                "end": {"value": 18, "modifier": "EXACT"},
            },
            "description": "Signal peptide",
            "evidences": [{"evidenceCode": "ECO:0000269"}],
        },
        {
            "type": "Disulfide bond",
            "location": {
                "start": {"value": 77, "modifier": "EXACT"},
                "end": {"value": 86, "modifier": "EXACT"},
            },
            "description": "Intrachain disulfide",
            "evidences": [{"evidenceCode": "ECO:0000269"}],
        },
        {
            "type": "Transmembrane",
            "location": {
                "start": {"value": 120, "modifier": "EXACT"},
                "end": {"value": 140, "modifier": "EXACT"},
            },
            "description": "Helical transmembrane",
            "evidences": [{"evidenceCode": "ECO:0000255"}],
        },
    ],
}


def test_uniprot_fetcher_with_mock() -> None:
    fetcher = UniProtSidecarFetcher()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P99999")

    assert sidecar["status"] == "AVAILABLE"
    assert sidecar["reviewed"] is False  # Mandatory human review invariant
    assert sidecar["primaryAccession"] == "P99999"
    assert sidecar["protein_length"] == 606
    assert sidecar["biological_validation"] == "NOT_AVAILABLE"
    assert len(sidecar["features"]) == 3
    assert len(sidecar["candidate_regions"]) == 3
    # Recommended dwell positions: signal cleavage (18), disulfide C2 (86), TM boundaries (120, 140)
    assert 18 in sidecar["recommended_positions"]
    assert 86 in sidecar["recommended_positions"]
    assert 120 in sidecar["recommended_positions"]
    assert 140 in sidecar["recommended_positions"]


def test_uniprot_sidecar_nbfold_contract_compliance() -> None:
    """Verify that generated sidecar passes nbfold contract validation without error."""
    fetcher = UniProtSidecarFetcher()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P99999")

    seq = sidecar["sequence"]["value"]
    seq_sha256 = sidecar["protein_sha256"]

    # Must pass nbfold.annotations.analyze_uniprot_tm_signal_annotations
    res = analyze_uniprot_tm_signal_annotations(
        protein_sequence=seq,
        protein_sha256=seq_sha256,
        uniprot_sidecar_data=sidecar,
        uniprot_sidecar_sha256=None,
        recommended_positions=[18, 86, 120, 140],
    )

    assert res["status"] in ("AVAILABLE", "CONFLICT")
    assert res["tool"] == "uniprot_sidecar_annotator"
    assert len(res["features"]) == 3
    assert res["accession"] == "P99999"


def test_uniprot_fetcher_sequence_mismatch_guard() -> None:
    fetcher = UniProtSidecarFetcher()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P99999", target_protein_sequence="ACDEF")

    assert sidecar["status"] == "ERROR_SEQUENCE_MISMATCH"
    assert "does not match" in sidecar["reason"]


def test_uniprot_fetcher_invalid_accession_guard() -> None:
    fetcher = UniProtSidecarFetcher()
    res1 = fetcher.fetch_and_format_sidecar("")
    assert res1["status"] == "ERROR_INVALID_ACCESSION"

    res2 = fetcher.fetch_and_format_sidecar("$$$BAD_ACC$$$")
    assert res2["status"] == "ERROR_INVALID_ACCESSION"


def test_uniprot_fetcher_empty_target_sequence_guard() -> None:
    fetcher = UniProtSidecarFetcher()
    res = fetcher.fetch_and_format_sidecar("P99999", target_protein_sequence="   ")
    assert res["status"] == "ERROR_INVALID_SEQUENCE"


def test_uniprot_fetcher_no_version_hallucination() -> None:
    """When entryAudit is omitted from raw API response, do NOT invent arbitrary version 1."""
    no_audit_record = dict(SAMPLE_MOCK_UNIPROT_RECORD)
    del no_audit_record["entryAudit"]

    fetcher = UniProtSidecarFetcher()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=no_audit_record,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P99999")

    assert sidecar["entryAudit"] is None


def test_uniprot_fetcher_out_of_bounds_features_safely_ignored() -> None:
    """Features with coordinates outside sequence bounds or start > end must be ignored."""
    corrupted_features_record = {
        "primaryAccession": "P99999",
        "sequence": {"value": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSE", "length": 30},
        "features": [
            {
                "type": "Signal",
                "location": {"start": {"value": 1}, "end": {"value": 15}},
                "description": "Valid signal",
            },
            {
                "type": "Transmembrane",
                "location": {"start": {"value": 25}, "end": {"value": 50}},  # End exceeds length 30
                "description": "Out of bounds",
            },
            {
                "type": "Disulfide bond",
                "location": {"start": {"value": 20}, "end": {"value": 10}},  # Start > End
                "description": "Inverted coordinates",
            },
        ],
    }

    fetcher = UniProtSidecarFetcher()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=corrupted_features_record,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P99999")

    assert len(sidecar["features"]) == 1
    assert sidecar["features"][0]["description"] == "Valid signal"


def test_uniprot_fetcher_cache_mismatch_safety(tmp_path: Path) -> None:
    """If a cache file exists with a different accession, it must not be returned."""
    cache_dir = tmp_path / "uniprot_cache"
    cache_dir.mkdir()
    bad_cache = cache_dir / "P99999_uniprot.json"
    bad_cache.write_text('{"primaryAccession": "P11111", "sequence": {"value": "ACD"}}', encoding="utf-8")

    fetcher = UniProtSidecarFetcher(cache_dir=cache_dir)
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,
    ):
        sidecar = fetcher.fetch_and_format_sidecar("P99999")

    assert sidecar["primaryAccession"] == "P99999"
    assert sidecar["protein_length"] == 606
