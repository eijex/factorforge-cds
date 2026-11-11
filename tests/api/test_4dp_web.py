"""Synthetic direct 4DP API tests with actual DPV240Optimizer and sidecar."""

import json
from io import BytesIO

from Bio.Seq import Seq
import pytest

from api.optimize import handler

META_PROTEIN = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"


def meta_sample_sidecar():
    return handler.generate_meta_sample_sidecar(None)


def post(data, *, path='/api/optimize', address='127.0.0.1', headers=None):
    h = object.__new__(handler)
    body = json.dumps(data).encode()
    h.path = path
    h.client_address = (address, 1234)
    h.headers = {"Content-Length": str(len(body)), "Host": "127.0.0.1", **(headers or {})}
    h.rfile = BytesIO(body)
    h.wfile = BytesIO()
    codes = []
    h.send_response = codes.append
    h.send_header = lambda *a: None
    h.end_headers = lambda: None
    h.do_POST()
    return codes[0], json.loads(h.wfile.getvalue())


@pytest.fixture
def local_preview(monkeypatch):
    monkeypatch.delenv("VERCEL", raising=False)
    monkeypatch.setenv("FACTORFORGE_LOCAL_4DP_FIXTURES", "true")


def test_4dp_sample_sidecar_action(local_preview):
    code, res = post({"action": "get_sample_sidecar"})
    assert code == 200
    assert res["success"] is True
    assert "sidecar" in res
    assert res["sidecar"]["pause_candidates"][0]["tunnel_pause_position"] == 38


def test_4dp_web_optimization_mock_blocking(local_preview):
    sidecar = meta_sample_sidecar()
    code, res = post({
        "sequence": META_PROTEIN,
        "design_mode": "paired_4dp",
        "pause_sidecar": sidecar,
        "positions_reviewed": True,
        "allow_mock": False,
    })
    # Should fail because mock sidecar is blocked without allow_mock=True
    assert code == 400
    assert res["success"] is False


def test_4dp_web_optimization_success(local_preview):
    sidecar = meta_sample_sidecar()
    code, res = post({
        "sequence": META_PROTEIN,
        "design_mode": "paired_4dp",
        "pause_sidecar": sidecar,
        "positions_reviewed": True,
        "allow_mock": True,
        "constraints": {"gc_min": 35, "gc_max": 55},
    })
    assert code == 200
    assert res["success"] is True
    assert res["design_mode"] == "paired_4dp"
    assert res["pause_positions"] == [38]
    assert len(res["candidates"]) == 2
    assert res["candidates"][0]["id"] == "baseline_1dp"
    assert res["candidates"][1]["id"] == "structure_4dp"
    assert res["recommended_candidate"]["id"] == "structure_4dp"
    assert res["metrics"]["cai"] > 0
    assert res["constraint_report"]["passed"] is True
    assert res["exit_tunnel_offset"] == 35
    assert res["model_protein_length"] == len(META_PROTEIN)
    assert all(str(Seq(c["dna_sequence"]).translate()) == META_PROTEIN for c in res["candidates"])
    assert "offset_evidence" in res
    if res["offset_evidence"]:
        assert "verified_facts" in res["offset_evidence"]
        assert "hypothesized_effects" in res["offset_evidence"]


@pytest.mark.parametrize("fields", [
    {"design_mode": "paired_4dp"},
    {"pause_sidecar": {}},
    {"sidecar_data": {}},
    {"action": "get_sample_sidecar"},
])
def test_public_preview_gate(monkeypatch, fields):
    monkeypatch.setattr("api.optimize.PUBLIC_4DP_ENABLED", False)
    status, result = post(fields)
    assert status == 403
    assert result["error_code"] == "EXPERIMENTAL_4DP_DISABLED"


def test_public_prepared_sidecar_runs_without_fixture_mode(monkeypatch):
    monkeypatch.delenv("FACTORFORGE_LOCAL_4DP_FIXTURES", raising=False)
    sidecar = meta_sample_sidecar()
    sidecar["is_mock"] = False
    sidecar["source_type"] = "ESMFOLD_PREPARED"
    code, res = post(
        {
            "sequence": META_PROTEIN,
            "design_mode": "paired_4dp",
            "pause_sidecar": sidecar,
            "positions_reviewed": True,
            "constraints": {"gc_min": 35, "gc_max": 55},
        },
        address="192.0.2.1",
        headers={"Host": "factorforge.eijex.com"},
    )
    assert code == 200
    assert res["success"] is True
    assert res["design_mode"] == "paired_4dp"
    assert res["pause_positions"] == [38]
    # Synthetic fixture endpoint remains forbidden without local fixture opt-in
    fixture_code, _ = post({"action": "get_sample_sidecar"}, headers={"Host": "127.0.0.1"})
    assert fixture_code == 403


def test_preview_gate_blocks_external_host_and_vercel(local_preview, monkeypatch):
    status, _ = post({"action": "get_sample_sidecar"}, address="192.0.2.1")
    assert status == 403
    status, _ = post({"action": "get_sample_sidecar"}, headers={"Origin": "https://example.test"})
    assert status == 403
    status, _ = post({"action": "get_sample_sidecar"}, headers={"Host": "example.test"})
    assert status == 403
    monkeypatch.setenv("VERCEL", "1")
    status, _ = post({"action": "get_sample_sidecar"})
    assert status == 403


def test_preview_rejects_non_inline_and_wrong_route(local_preview):
    fields = {
        "sequence": META_PROTEIN, "design_mode": "paired_4dp",
        "pause_sidecar": "must-not-be-opened.json", "positions_reviewed": True,
    }
    status, result = post(fields)
    assert status == 400
    assert "inline JSON object" in result["error"]
    status, _ = post(fields, path="/api/optimize/batch")
    assert status == 400


@pytest.mark.parametrize("changes", [
    {"positions_reviewed": False},
    {"exit_tunnel_offset": 42},
    {"sequence": "A" + META_PROTEIN[1:]},
    {"sequence": "M" * 501},
    {"allow_mock": "true"},
    {"constraints": {"gc_min": 90, "gc_max": 10}},
])
def test_preview_rejects_invalid_or_unreviewed_input(local_preview, changes):
    fields = {
        "sequence": META_PROTEIN, "design_mode": "paired_4dp",
        "pause_sidecar": meta_sample_sidecar(), "positions_reviewed": True, "allow_mock": True,
    }
    status, result = post({**fields, **changes})
    assert status == 400
    assert result["success"] is False


def test_4dp_web_accepts_472aa_synthetic_sequence(monkeypatch):
    from factorforge.structure.alphafold_server_ingest import _pause_candidates

    monkeypatch.delenv("FACTORFORGE_LOCAL_4DP_FIXTURES", raising=False)
    protein_472 = "M" + ("ACDEFGHIKLMNPQRSTVWY" * 23) + "GALMQRSTVWY"
    residues = [
        {"seq_id": i, "aa": aa, "plddt": 45.0 if 100 <= i <= 115 else 92.0}
        for i, aa in enumerate(protein_472, 1)
    ]
    pause_candidates = _pause_candidates(residues, 35)
    assert pause_candidates
    sidecar = {
        "status": "AVAILABLE",
        "source_type": "ESMFOLD_PREPARED",
        "is_mock": False,
        "reviewed": False,
        "coordinate_system": "PDB_res_seq_1based",
        "primary_chain": "A",
        "tunnel_offset_applied": 35,
        "pause_candidates": pause_candidates,
        "chains": {
            "A": {
                "chain_id": "A",
                "sequence": protein_472,
                "length": 472,
                "residues": residues,
                "pause_candidates": pause_candidates,
            }
        },
    }
    code, res = post(
        {
            "sequence": protein_472,
            "design_mode": "paired_4dp",
            "pause_sidecar": sidecar,
            "positions_reviewed": True,
            "constraints": {"gc_min": 35, "gc_max": 55},
        },
        address="192.0.2.1",
        headers={"Host": "factorforge.eijex.com"},
    )
    assert code == 200, res
    assert res["model_protein_length"] == 472
    assert len(res["recommended_candidate"]["dna_sequence"]) == 472 * 3
