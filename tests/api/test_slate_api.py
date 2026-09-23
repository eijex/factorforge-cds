"""Integration tests for POST /api/slate API handler (Job 283A)."""

import json
from api.optimize import handler


class MockWFile:
    def __init__(self):
        self.data = b""

    def write(self, b):
        self.data += b


def call_api(path: str, body: dict, method: str = "POST"):
    h = handler.__new__(handler)
    h.path = path
    h.command = method
    h.wfile = MockWFile()
    h.headers = {"Content-Type": "application/json"}
    h._headers_buffer = []

    def mock_send_response(code, message=None):
        h.status_code = code

    def mock_send_header(k, v):
        h._headers_buffer.append((k, v))

    def mock_end_headers():
        pass

    h.send_response = mock_send_response
    h.send_header = mock_send_header
    h.end_headers = mock_end_headers

    if method == "GET":
        h.do_GET()
    else:
        # Mock reading request body
        import io

        body_bytes = json.dumps(body).encode("utf-8")
        h.rfile = io.BytesIO(body_bytes)
        h.headers["Content-Length"] = str(len(body_bytes))
        h.do_POST()

    response_json = json.loads(h.wfile.data.decode("utf-8"))
    return h.status_code, response_json


def test_api_slate_endpoint_success():
    payload = {
        "sequence": "DIQMTQSPSSLSASVGDRVT",
        "target_name": "Humira-LC-Fragment",
        "top_k": 3,
        "novelty_class": "Class_A_InDistribution",
    }
    status, res = call_api("/api/slate", payload)
    assert status == 200
    assert res["success"] is True
    data = res["data"]
    assert data["$schema"] == "https://eijex.com/schemas/factorforge/candidate-slate-v0.1.json"
    assert data["target_metadata"]["target_name"] == "Humira-LC-Fragment"
    assert data["slate_summary"]["top_k_count"] <= 3
    assert len(data["candidates"]) >= 1


def test_api_optimize_endpoint_backward_compatibility():
    """Verify that /api/optimize remains 100% semantically unchanged for single sequence requests."""
    payload = {
        "sequence": "DIQMTQSPSSLSASVGDRVT",
        "profile": "balanced",
    }
    status, res = call_api("/api/optimize", payload)
    assert status == 200
    assert "optimized_sequence" in res
    assert "metrics" in res
    assert "cai" in res["metrics"]
    assert "gc_percent" in res["metrics"]


def test_api_slate_v2_endpoint_success():
    """Verify POST /api/slate Slate v2 (Job 293) request returns exact K=25 slate with validation summary."""
    payload = {
        "protein_sequence": "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAK",
        "host": "nbenthamiana",
        "target_gc": 0.45,
        "slate_size": 25,
        "stop_policy": "append_preferred",
        "generation_mode": "deterministic_beam",
        "weights": {
            "cai": 0.40,
            "gc_fidelity": 0.25,
            "mfe_initiation": 0.20,
            "rare_codon_guard": 0.15,
        },
        "forbidden_sites": {
            "type_iis": ["BsaI", "BsmBI"],
            "additional": ["NotI", "XhoI"],
        },
    }
    status, res = call_api("/api/slate", payload)
    assert status == 200
    assert res["status"] == "success"
    assert "job_id" in res
    assert res["metadata"]["slate_size"] == 25
    assert len(res["slate"]) == 25
    assert res["validation_summary"]["aa_conservation_by_construction"] is True
    assert res["validation_summary"]["aa_conservation_verified"] is True
    assert res["validation_summary"]["aa_identity_pct"] == 100.0
    assert res["validation_summary"]["forbidden_sites_clean"] is True
    assert res["validation_summary"]["primary_host_status"] == "PASS"

    first_cand = res["slate"][0]
    assert first_cand["rank"] == 1
    assert "cross_host_sensitivity" in first_cand
    assert first_cand["cross_host_sensitivity"]["primary_host"]["cai_reference"] == 1.00


def test_api_slate_v2_fail_closed_validation():
    """Verify POST /api/slate returns HTTP 400 for invalid parameter contracts."""
    # 1. Invalid slate_size
    status, res = call_api("/api/slate", {
        "protein_sequence": "MKWVTFISLLLLFSSAYSRG",
        "slate_size": 999,
    })
    assert status == 400
    assert res["status"] == "error"

    # 2. Invalid stop_policy
    status, res = call_api("/api/slate", {
        "protein_sequence": "MKWVTFISLLLLFSSAYSRG",
        "stop_policy": "invalid_policy_123",
    })
    assert status == 400
    assert res["status"] == "error"
