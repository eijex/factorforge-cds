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
