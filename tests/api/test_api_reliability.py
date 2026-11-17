"""Reliability tests for /api/optimize and scripts/serve.py (Job 681 Phase 2)."""

from __future__ import annotations

from io import BytesIO
import json
from unittest.mock import patch

import pytest

import api.optimize as optimize_api
from api.optimize import handler


@pytest.fixture(autouse=True)
def _clean_api_state() -> None:
    optimize_api._reset_rate_limiter()
    optimize_api._reset_esmfold_cache()
    yield
    optimize_api._reset_rate_limiter()
    optimize_api._reset_esmfold_cache()


def _post_raw(
    path: str,
    raw_body: bytes,
    headers: dict[str, str] | None = None,
    client_ip: str | None = "203.0.113.10",
) -> tuple[int, dict[str, str], dict]:
    request_handler = object.__new__(handler)
    request_handler.path = path
    request_headers = {"Content-Length": str(len(raw_body))}
    if headers is not None:
        request_headers = dict(headers)
    request_handler.headers = request_headers
    if client_ip is not None:
        request_handler.client_address = (client_ip, 54321)
    request_handler.rfile = BytesIO(raw_body)
    request_handler.wfile = BytesIO()
    responses: list[int] = []
    out_headers: dict[str, str] = {}
    request_handler.send_response = responses.append
    request_handler.send_header = lambda key, value: out_headers.__setitem__(key, value)
    request_handler.end_headers = lambda: None

    handler.do_POST(request_handler)

    request_handler.wfile.seek(0)
    response_body = json.loads(request_handler.wfile.read().decode("utf-8"))
    return responses[0], out_headers, response_body


def _post_json(
    path: str,
    payload: dict,
    *,
    headers: dict[str, str] | None = None,
    client_ip: str | None = "203.0.113.10",
) -> tuple[int, dict[str, str], dict]:
    body = json.dumps(payload).encode("utf-8")
    merged_headers = {"Content-Length": str(len(body))}
    if headers:
        merged_headers.update(headers)
    return _post_raw(path, body, headers=merged_headers, client_ip=client_ip)


def test_rate_limit_returns_429_with_retry_after(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(optimize_api, "RATE_LIMIT_MAX_REQUESTS", 2)
    monkeypatch.setattr(optimize_api, "RATE_LIMIT_WINDOW_SECONDS", 60)

    payload = {"sequence": "MSKGEELFTGVVPILVELD", "profile": "balanced"}
    s1, _, _ = _post_json("/api/optimize", payload, client_ip="198.51.100.1")
    s2, _, _ = _post_json("/api/optimize", payload, client_ip="198.51.100.1")
    s3, hdrs3, body3 = _post_json("/api/optimize", payload, client_ip="198.51.100.1")

    assert s1 == 200
    assert s2 == 200
    assert s3 == 429
    assert body3["success"] is False
    assert body3["error_code"] == "RATE_LIMIT_EXCEEDED"
    assert int(body3["retry_after"]) >= 1
    assert int(hdrs3["Retry-After"]) >= 1

    # Independent IP bucket is unaffected
    s_other, _, _ = _post_json("/api/optimize", payload, client_ip="198.51.100.2")
    assert s_other == 200


def test_rate_limit_respects_x_forwarded_for_and_window_expiry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(optimize_api, "RATE_LIMIT_MAX_REQUESTS", 1)
    monkeypatch.setattr(optimize_api, "RATE_LIMIT_WINDOW_SECONDS", 10)

    allowed1, _ = optimize_api._check_rate_limit("192.0.2.55", now=100.0)
    allowed2, retry_after = optimize_api._check_rate_limit("192.0.2.55", now=102.0)
    allowed3, _ = optimize_api._check_rate_limit("192.0.2.55", now=111.0)

    assert allowed1 is True
    assert allowed2 is False
    assert retry_after >= 1
    assert allowed3 is True

    # Verify X-Forwarded-For takes precedence over socket client_address
    payload = {"sequence": "MSKGEELFTGVVPILVELD", "profile": "balanced"}
    s_xff1, _, _ = _post_json(
        "/api/optimize",
        payload,
        headers={"X-Forwarded-For": "203.0.113.99, 10.0.0.1"},
        client_ip="127.0.0.1",
    )
    s_xff2, _, _ = _post_json(
        "/api/optimize",
        payload,
        headers={"X-Forwarded-For": "203.0.113.99, 10.0.0.2"},
        client_ip="127.0.0.1",
    )
    assert s_xff1 == 200
    assert s_xff2 == 429


def test_esmfold_sha256_lru_cache_avoids_duplicate_remote_calls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(optimize_api, "ESMFOLD_CACHE_MAX_ENTRIES", 2)
    fake_resolved = {
        "content": "HEADER FAKE PDB\n",
        "format": "pdb",
        "chain": "A",
        "model": 1,
        "positions": [12],
        "positions_reviewed": False,
        "source_label": "ESMFold (api.esmatlas.com · UNVERIFIED_PUBLIC_CORPUS)",
        "protein_sha256": "abc",
        "structure_sha256": "def",
        "ca_trace": [],
        "cartoon_pdb": "",
        "sidecar_4dp": None,
    }

    with patch(
        "factorforge.structure.auto_resolver.auto_resolve_structure_for_protein",
        return_value=fake_resolved,
    ) as mock_resolve:
        s1, _, b1 = _post_json(
            "/api/optimize",
            {
                "action": "auto_resolve_structure",
                "sequence": "MSKGEELFTGVVPILVELD",
                "allow_remote": True,
            },
        )
        # Same protein with lowercase/whitespace/stop should hit the same SHA-256 entry
        s2, _, b2 = _post_json(
            "/api/optimize",
            {
                "action": "auto_resolve_structure",
                "sequence": " mskgeelftgvvpilveld* ",
                "allow_remote": True,
            },
        )

    assert s1 == 200
    assert b1["cache_hit"] is False
    assert s2 == 200
    assert b2["cache_hit"] is True
    assert b2["resolved"]["content"] == "HEADER FAKE PDB\n"
    mock_resolve.assert_called_once()

    # Verify LRU eviction when exceeding ESMFOLD_CACHE_MAX_ENTRIES=2
    optimize_api._set_cached_esmfold("MAAAA", {"id": "a"})
    optimize_api._set_cached_esmfold("MCCCC", {"id": "c"})
    # Touch MAAAA so MCCCC becomes LRU
    assert optimize_api._get_cached_esmfold("MAAAA") == {"id": "a"}
    optimize_api._set_cached_esmfold("MDDDD", {"id": "d"})
    assert optimize_api._get_cached_esmfold("MCCCC") is None
    assert optimize_api._get_cached_esmfold("MAAAA") == {"id": "a"}
    assert optimize_api._get_cached_esmfold("MDDDD") == {"id": "d"}


def test_batch_and_compare_enforce_total_aa_caps(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(optimize_api, "MAX_BATCH_TOTAL_AA", 50)
    monkeypatch.setattr(optimize_api, "MAX_COMPARE_TOTAL_AA", 50)

    seq_30 = "M" + "A" * 29
    batch_status, _, batch_body = _post_json(
        "/api/optimize/batch",
        {
            "sequences": [
                {"id": "s1", "sequence": seq_30},
                {"id": "s2", "sequence": seq_30},
            ],
            "profile": "balanced",
        },
    )
    assert batch_status == 400
    assert batch_body["success"] is False
    assert "amino-acid equivalents" in batch_body["error"]

    compare_status, _, compare_body = _post_json(
        "/api/optimize/compare",
        {
            "sequence": seq_30,
            "profiles": ["balanced", "high_cai"],
        },
    )
    assert compare_status == 400
    assert compare_body["success"] is False
    assert "total amino-acid equivalents" in compare_body["error"]


@pytest.mark.parametrize(
    ("headers", "raw_body", "expected_status", "expected_substr"),
    [
        ({"Content-Length": "abc"}, b"{}", 400, "Invalid Content-Length header"),
        ({"Content-Length": "-5"}, b"{}", 400, "Invalid Content-Length header"),
        ({"Content-Length": "2000000"}, b"{}", 413, "1 MB limit"),
        ({"Transfer-Encoding": "chunked", "Content-Length": "2"}, b"{}", 400, "Chunked"),
        ({"Content-Length": "2"}, b"\xff\xfe", 400, "valid UTF-8"),
    ],
)
def test_do_post_hardened_content_length_and_body_parsing(
    headers: dict[str, str],
    raw_body: bytes,
    expected_status: int,
    expected_substr: str,
) -> None:
    status, _, body = _post_raw("/api/optimize", raw_body, headers=headers)
    assert status == expected_status
    assert body["success"] is False
    assert expected_substr in body["error"]


def test_scripts_serve_uses_threading_http_server() -> None:
    import importlib.util
    from http.server import ThreadingHTTPServer
    from pathlib import Path

    serve_path = Path(__file__).resolve().parents[2] / "scripts" / "serve.py"
    spec = importlib.util.spec_from_file_location("factorforge_serve", serve_path)
    assert spec is not None and spec.loader is not None
    source = serve_path.read_text(encoding="utf-8")
    assert "ThreadingHTTPServer" in source
    assert "HTTPServer(" not in source.replace("ThreadingHTTPServer(", "")
    assert ThreadingHTTPServer is not None
