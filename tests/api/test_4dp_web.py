"""Synthetic direct 4DP API tests with actual DPV240Optimizer and sidecar."""

import json
from io import BytesIO

from api.optimize import handler

META_PROTEIN = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"


def meta_sample_sidecar():
    return handler.generate_meta_sample_sidecar(None)


def post(data):
    h = object.__new__(handler)
    body = json.dumps(data).encode()
    h.path = "/api/optimize"
    h.headers = {"Content-Length": str(len(body))}
    h.rfile = BytesIO(body)
    h.wfile = BytesIO()
    codes = []
    h.send_response = codes.append
    h.send_header = lambda *a: None
    h.end_headers = lambda: None
    h.do_POST()
    return codes[0], json.loads(h.wfile.getvalue())


def test_4dp_sample_sidecar_action():
    code, res = post({"action": "get_sample_sidecar"})
    assert code == 200
    assert res["success"] is True
    assert "sidecar" in res
    assert res["sidecar"]["pause_candidates"][0]["tunnel_pause_position"] == 38


def test_4dp_web_optimization_mock_blocking():
    sidecar = meta_sample_sidecar()
    code, res = post({
        "sequence": META_PROTEIN,
        "design_mode": "paired_4dp",
        "sidecar_data": sidecar,
        "allow_mock": False,
    })
    # Should fail because mock sidecar is blocked without allow_mock=True
    assert code in (400, 500) or res.get("success") is False


def test_4dp_web_optimization_success():
    sidecar = meta_sample_sidecar()
    code, res = post({
        "sequence": META_PROTEIN,
        "design_mode": "paired_4dp",
        "sidecar_data": sidecar,
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
    assert "sllm_report" in res
    if res["sllm_report"]:
        assert "verified_facts" in res["sllm_report"]
        assert "hypothesized_effects" in res["sllm_report"]
