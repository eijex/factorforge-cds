"""Unit tests for ESMFold client, PDB parsing, and sidecar generation."""

import pytest

from factorforge.structure.esmfold_client import (
    ESMFoldClientError,
    _clean_sequence,
    _generate_mock_pdb,
    format_esm_report,
    parse_esm_pdb,
    predict_and_ingest_esm,
    predict_structure_esm,
)

SAMPLE_SEQ = "MKIFLMMMKIFLMMMKIFLMMMKIFLMMMKIFLMM"


def test_clean_sequence_valid():
    assert _clean_sequence("  mkif \n lmm* ") == "MKIFLMM"


def test_clean_sequence_empty():
    with pytest.raises(ESMFoldClientError, match="cannot be empty"):
        _clean_sequence("   *  ")


def test_clean_sequence_invalid_aa():
    with pytest.raises(ESMFoldClientError, match="non-standard amino acid"):
        _clean_sequence("MKIFBZX")


def test_generate_mock_pdb():
    pdb = _generate_mock_pdb("MKIFLMM")
    lines = pdb.splitlines()
    assert lines[0].startswith("HEADER")
    ca_lines = [line for line in lines if line.startswith("ATOM") and " CA " in line]
    assert len(ca_lines) == 7
    assert lines[-1] == "END"


def test_parse_esm_pdb_mock():
    pdb = _generate_mock_pdb(SAMPLE_SEQ)
    result = parse_esm_pdb(pdb, target_chain="A", tunnel_offset=35)

    assert result["status"] == "AVAILABLE"
    assert result["source_type"] == "ESMFold"
    assert result["coordinate_system"] == "PDB_res_seq_1based"
    assert result["primary_chain"] == "A"
    assert result["tunnel_offset_applied"] == 35

    chain = result["chains"]["A"]
    assert chain["length"] == len(SAMPLE_SEQ)
    assert chain["sequence"] == SAMPLE_SEQ
    assert len(chain["residues"]) == len(SAMPLE_SEQ)

    # First residue should be 1-based Met
    assert chain["residues"][0]["seq_id"] == 1
    assert chain["residues"][0]["aa"] == "M"

    # There should be pause candidates (e.g., residues 1..3 have plddt 45.0, so end 3 + 35 = 38)
    candidates = chain["pause_candidates"]
    assert len(candidates) >= 1
    assert candidates[0]["start_residue"] == 1
    assert candidates[0]["end_residue"] == 3
    assert candidates[0]["tunnel_pause_position"] == 38


def test_parse_esm_pdb_empty():
    with pytest.raises(ESMFoldClientError, match="content is empty"):
        parse_esm_pdb("")


def test_parse_esm_pdb_missing_ca():
    pdb = "ATOM      1  N   MET A   1      1.0 1.0 1.0  1.00 80.00           N\nEND\n"
    with pytest.raises(ESMFoldClientError, match="No CA atoms found"):
        parse_esm_pdb(pdb)


def test_parse_esm_pdb_sequence_mismatch():
    pdb = _generate_mock_pdb("MKIFLMM")
    with pytest.raises(ESMFoldClientError, match="does not match expected sequence"):
        parse_esm_pdb(pdb, expected_sequence="MKIFLMA")



def test_predict_structure_esm_mock():
    pdb = predict_structure_esm(SAMPLE_SEQ, provider="mock")
    assert "HEADER    MOCK ESMFOLD PREDICTION" in pdb
    assert " CA " in pdb


def test_predict_structure_esm_unknown_provider():
    with pytest.raises(ESMFoldClientError, match="Unknown provider"):
        predict_structure_esm(SAMPLE_SEQ, provider="invalid_provider")


def test_predict_and_ingest_esm(tmp_path):
    sidecar_path = tmp_path / "sidecar.json"
    pdb_path = tmp_path / "structure.pdb"

    result, pdb_text = predict_and_ingest_esm(
        sequence=SAMPLE_SEQ,
        provider="mock",
        output_sidecar_path=sidecar_path,
        output_pdb_path=pdb_path,
        tunnel_offset=35,
    )

    assert sidecar_path.is_file()
    assert pdb_path.is_file()
    assert result["status"] == "AVAILABLE"
    assert result["chains"]["A"]["length"] == len(SAMPLE_SEQ)
    assert pdb_path.read_text(encoding="utf-8") == pdb_text


def test_format_esm_report():
    pdb = _generate_mock_pdb(SAMPLE_SEQ)
    result = parse_esm_pdb(pdb, tunnel_offset=35)
    report = format_esm_report(result)
    assert "ESMFold Structure Report" in report
    assert f"Protein Length: {len(SAMPLE_SEQ)} aa" in report
    assert "Tunnel Pause Position: 38" in report


def test_mismatched_model_cannot_write_artifacts(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "factorforge.structure.esmfold_client.predict_structure_esm",
        lambda **kwargs: _generate_mock_pdb("LKIF"),
    )
    sidecar = tmp_path / "sidecar.json"
    pdb = tmp_path / "model.pdb"
    with pytest.raises(ESMFoldClientError, match="does not match"):
        predict_and_ingest_esm("MKIF", output_sidecar_path=sidecar, output_pdb_path=pdb)
    assert not sidecar.exists()
    assert not pdb.exists()


def test_existing_outputs_are_preserved(tmp_path):
    output = tmp_path / "existing.pdb"
    output.write_text("original", encoding="utf-8")
    sidecar = tmp_path / "sidecar.json"
    with pytest.raises(ESMFoldClientError, match="Output already exists"):
        predict_and_ingest_esm(
            "MKIF", provider="mock", output_pdb_path=output, output_sidecar_path=sidecar
        )
    assert output.read_text(encoding="utf-8") == "original"
    assert not sidecar.exists()


@pytest.mark.parametrize("value", ["   nan", "   inf", "101.00", " -1.00"])
def test_invalid_confidence_is_rejected(value):
    lines = _generate_mock_pdb("MKIF").splitlines()
    lines[2] = lines[2][:60] + value + lines[2][66:]
    with pytest.raises(ESMFoldClientError, match="CA confidence"):
        parse_esm_pdb("\n".join(lines))


def test_hosted_nvidia_requires_configured_endpoint(monkeypatch):
    monkeypatch.delenv("NVIDIA_ESMFOLD_URL", raising=False)
    with pytest.raises(ESMFoldClientError, match="NVIDIA_ESMFOLD_URL"):
        predict_structure_esm("MKIF", provider="nvidia")


def test_remote_contract_disables_redirects(monkeypatch):
    from types import SimpleNamespace

    def post(url, **kwargs):
        assert url == "https://example.test/esmfold"
        assert kwargs["data"] == "MKIF"
        assert kwargs["allow_redirects"] is False
        return SimpleNamespace(status_code=302)

    monkeypatch.setenv("ESM_ATLAS_URL", "https://example.test/esmfold")
    monkeypatch.setattr("factorforge.structure.esmfold_client.requests.post", post)
    with pytest.raises(ESMFoldClientError, match="HTTP 302"):
        predict_structure_esm("MKIF")
