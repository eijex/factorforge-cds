from __future__ import annotations

import json
import sys
from io import BytesIO
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from api.optimize import handler


def _call_optimize(payload: dict) -> tuple[int, dict]:
    body = json.dumps(payload).encode("utf-8")

    class _DummyHandler(handler):
        def __init__(self):
            self.path = "/api/optimize"
            self.headers = {"Content-Length": str(len(body)), "Host": "localhost:3000"}
            self.rfile = BytesIO(body)
            self.wfile = BytesIO()
            self.status_code = 200

        def send_response(self, code, message=None):
            self.status_code = code

        def send_header(self, keyword, value):
            pass

        def end_headers(self):
            pass

    h = _DummyHandler()
    h.do_POST()
    return h.status_code, json.loads(h.wfile.getvalue().decode("utf-8"))


def test_experimental_2dp_web_endpoint() -> None:
    status, res = _call_optimize(
        {
            "design_mode": "experimental_2dp",
            "sequence": "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGGSLRLSCAASGFTFD",
        }
    )
    assert status == 200
    assert res["success"] is True
    assert res["design_mode"] == "experimental_2dp"
    assert res["biological_validation"] == "NOT_AVAILABLE"
    assert "cpg_count" in res["metrics"]
    assert "tpa_count" in res["metrics"]


def test_public_5dp_request_is_rejected() -> None:
    status, _ = _call_optimize({"design_mode": "experimental_5dp"})
    assert status == 403
