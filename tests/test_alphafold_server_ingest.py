"""Unit tests for FactorForge AlphaFold Server structure ingest bridge."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from factorforge.cli.main import cli
from factorforge.structure.alphafold_server_ingest import (
    AlphaFoldIngestError,
    format_ingest_report,
    ingest_alphafold_server_result,
    parse_af3_cif,
    parse_af3_summary_confidences,
)

SAMPLE_CIF_PATH = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "factorforge"
    / "structure"
    / "cd9_alphafold3_model_0.cif"
)


def test_parse_af3_cif_live_sample():
    """Verify parsing live AlphaFold 3 mmCIF model with real pLDDT confidence."""
    assert SAMPLE_CIF_PATH.exists()
    res = parse_af3_cif(SAMPLE_CIF_PATH, target_chain="A")

    assert res["status"] == "AVAILABLE"
    assert res["chain_count"] == 1
    assert "A" in res["chains"]
    chain_a = res["chains"]["A"]

    assert chain_a["length"] == 236
    assert chain_a["mean_plddt"] > 80.0
    assert len(chain_a["sequence"]) == 236

    # Verify distribution sums to ~100%
    dist = chain_a["plddt_distribution"]
    total_pct = (
        dist["very_high_ge_90"]["percent"]
        + dist["confident_70_89"]["percent"]
        + dist["low_50_69"]["percent"]
        + dist["very_low_lt_50"]["percent"]
    )
    assert 99.0 <= total_pct <= 101.0

    # Verify pause candidate coordinates
    candidates = chain_a["pause_candidates"]
    assert len(candidates) > 0
    for c in candidates:
        assert c["tunnel_pause_position"] == c["end_residue"] + 35
        assert c["mean_plddt"] < 70.0


def test_parse_af3_summary_confidences(tmp_path: Path):
    """Verify extraction of global AF3 summary confidence metrics."""
    summary_file = tmp_path / "summary_confidences_0.json"
    summary_file.write_text(
        json.dumps({
            "ranking_score": 0.88,
            "ptm": 0.85,
            "iptm": 0.91,
            "fraction_disordered": 0.05,
            "has_clash": 0.0,
        }),
        encoding="utf-8",
    )

    metrics = parse_af3_summary_confidences(summary_file)
    assert metrics["ranking_score"] == 0.88
    assert metrics["ptm"] == 0.85
    assert metrics["iptm"] == 0.91
    assert metrics["fraction_disordered"] == 0.05
    assert metrics["has_clash"] is False


def test_ingest_alphafold_server_result_full_pipeline(tmp_path: Path):
    """Verify full ingestion pipeline writing structured sidecar."""
    summary_file = tmp_path / "summary.json"
    summary_file.write_text(
        json.dumps({"ranking_score": 0.82, "ptm": 0.79, "iptm": 0.84}),
        encoding="utf-8",
    )
    out_file = tmp_path / "af3_sidecar.json"

    res = ingest_alphafold_server_result(
        cif_path=SAMPLE_CIF_PATH,
        chain="A",
        summary_json_path=summary_file,
        output_path=out_file,
        tunnel_offset=35,
    )

    assert res["status"] == "AVAILABLE"
    assert out_file.exists()
    saved_data = json.loads(out_file.read_text(encoding="utf-8"))
    assert saved_data["status"] == "AVAILABLE"
    assert saved_data["summary_confidences"]["ptm"] == 0.79
    assert "epistemic_status" in saved_data

    # Verify format_ingest_report
    text_rep = format_ingest_report(res)
    assert "AlphaFold Server (AF3) Ingest Report" in text_rep
    assert "Chain 'A'" in text_rep


def test_cli_ingest_af3_text_and_json(tmp_path: Path):
    """Verify CLI structure ingest-af3 execution for text and JSON outputs."""
    runner = CliRunner()

    # Plain text run
    res_text = runner.invoke(
        cli,
        ["structure", "ingest-af3", "-c", str(SAMPLE_CIF_PATH), "--chain", "A"],
    )
    assert res_text.exit_code == 0
    assert "AlphaFold Server (AF3) Ingest Report" in res_text.output
    assert "Mean pLDDT:" in res_text.output
    assert "Epistemic: UNVERIFIED_PUBLIC_CORPUS" in res_text.output

    # JSON run
    res_json = runner.invoke(
        cli,
        ["structure", "ingest-af3", "-c", str(SAMPLE_CIF_PATH), "--json"],
    )
    assert res_json.exit_code == 0
    parsed = json.loads(res_json.output)
    assert parsed["status"] == "AVAILABLE"
    assert "chains" in parsed


def test_ingest_af3_error_handling(tmp_path: Path):
    """Verify descriptive error handling on corrupt or missing mmCIF data."""
    bad_cif = tmp_path / "corrupt.cif"
    bad_cif.write_text("data_test\n_entry.id test\n# Empty without atom_site", encoding="utf-8")

    with pytest.raises(AlphaFoldIngestError, match="No _atom_site loop found"):
        parse_af3_cif(bad_cif)

    runner = CliRunner()
    res_fail = runner.invoke(cli, ["structure", "ingest-af3", "-c", str(bad_cif)])
    assert res_fail.exit_code != 0
    assert "ingest error" in res_fail.output


def test_model_coordinates_do_not_bridge_missing_residues_or_certify_precursor(tmp_path):
    cif = tmp_path / "model.cif"
    cif.write_text(
        "data_test\nloop_\n_atom_site.label_atom_id\n_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n_atom_site.auth_seq_id\n"
        "_atom_site.label_comp_id\n_atom_site.B_iso_or_equiv\n"
        "'CA' A 1 101 ALA 42\nCA A 2 102 ALA 42\n"
        "CA A 5 105 ALA 42\nCA A 6 106 ALA 42\nCA A 7 107 ALA 42\n"
        "CA A 9 109 ALA 99\n", encoding="utf-8",
    )
    result = parse_af3_cif(cif, tunnel_offset=2)
    candidate, = result["chains"]["A"]["pause_candidates"]
    assert (candidate["start_residue"], candidate["end_residue"]) == (5, 7)
    assert candidate["tunnel_pause_position"] == 9
    assert candidate["within_structure_chain"]
    assert candidate["within_precursor"] is None

    original = cif.read_bytes()
    result = CliRunner().invoke(cli, ["structure", "ingest-af3", "-c", str(cif), "-o", str(cif)])
    assert result.exit_code != 0
    assert cif.read_bytes() == original


def test_invalid_summary_shapes_and_nonfinite_confidences():
    for summary in ('[]', '{"ptm": NaN}', '{"iptm": 1.1}'):
        with pytest.raises(AlphaFoldIngestError):
            parse_af3_summary_confidences(summary)
