"""Unit tests for FactorForge structure CLI commands."""

from pathlib import Path
from unittest.mock import patch
from click.testing import CliRunner
import json
from factorforge.cli.main import cli

SAMPLE_MOCK_UNIPROT_RECORD = {
    "primaryAccession": "P99999",
    "entryAudit": {"entryVersion": 10, "sequenceVersion": 2},
    "sequence": {
        "value": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSE",
        "length": 30,
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
    ],
}

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
                    "accession": "P21926",
                    "protein_length": 30,
                    "entry_protein_locations": [
                        {
                            "fragments": [
                                {
                                    "start": 5,
                                    "end": 25,
                                }
                            ]
                        }
                    ],
                }
            ],
        },
    ]
}

SAMPLE_PROTEIN_30 = "M" * 30


def test_structure_cli_fetch_uniprot() -> None:
    runner = CliRunner()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,
    ):
        result = runner.invoke(cli, ["structure", "fetch-uniprot", "-a", "P99999"])

    assert result.exit_code == 0
    assert "P99999" in result.output
    assert "CURATED_UNIPROT_RECORD" in result.output


def test_structure_cli_fetch_interpro() -> None:
    runner = CliRunner()
    with patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value=SAMPLE_MOCK_INTERPRO_RECORD,
    ):
        result = runner.invoke(cli, ["structure", "fetch-interpro", "-a", "P21926", "-s", SAMPLE_PROTEIN_30])

    assert result.exit_code == 0
    assert "IPR000452" in result.output
    assert "INTERPRO_REST_API" in result.output


def test_structure_cli_annotate_no_structure() -> None:
    runner = CliRunner()
    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,
    ), patch(
        "factorforge.structure.interpro_sidecar_fetcher.fetch_interpro_raw_record",
        return_value=SAMPLE_MOCK_INTERPRO_RECORD,
    ):
        result = runner.invoke(cli, ["structure", "annotate", "-a", "P99999", "--no-interpro"])

    assert result.exit_code == 0
    assert "viewer_annotations" in result.output
    assert "UNVERIFIED_PUBLIC_CORPUS" in result.output


def test_structure_cli_annotate_with_mismatched_structure() -> None:
    """When a structure file with mismatched sequence is provided, coordinate_mapping_verified must be False."""
    runner = CliRunner()
    cif_path = Path(__file__).resolve().parents[1] / "src" / "factorforge" / "structure" / "cd9_alphafold3_model_0.cif"

    with patch(
        "factorforge.structure.uniprot_sidecar_fetcher.fetch_uniprot_raw_record",
        return_value=SAMPLE_MOCK_UNIPROT_RECORD,  # P99999 sequence is 30 aa of MKWVT... (does not match CD9)
    ):
        result = runner.invoke(
            cli,
            ["structure", "annotate", "-a", "P99999", "--cif", str(cif_path), "--no-interpro"],
        )

    assert result.exit_code == 0
    assert "Warning: Structure coordinate mapping failed" in result.output
    # Parse json output from output string (take last JSON block)
    json_start = result.output.find("{")
    assert json_start != -1
    parsed = json.loads(result.output[json_start:])
    assert parsed["coordinate_mapping_verified"] is False
    assert len(parsed["viewer_annotations"]) > 0
    assert parsed["viewer_annotations"][0]["coordinate_mapping"] == "ASSUMED_1_TO_1_UNVERIFIED"
    assert parsed["viewer_annotations"][0]["mapping_status"] == "UNVERIFIED"
