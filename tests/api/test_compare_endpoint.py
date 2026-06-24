"""Tests for the /api/optimize/compare endpoint."""

from __future__ import annotations

from io import BytesIO
import json

from api.optimize import handler


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"


def _handler() -> handler:
    return object.__new__(handler)


def _post_compare(payload: dict) -> tuple[int, dict]:
    request_handler = _handler()
    body = json.dumps(payload).encode("utf-8")
    request_handler.path = "/api/optimize/compare"
    request_handler.headers = {"Content-Length": str(len(body))}
    request_handler.rfile = BytesIO(body)
    request_handler.wfile = BytesIO()
    responses: list[int] = []
    headers: list[tuple[str, str]] = []
    request_handler.send_response = responses.append
    request_handler.send_header = lambda key, value: headers.append((key, value))
    request_handler.end_headers = lambda: None

    handler.do_POST(request_handler)

    request_handler.wfile.seek(0)
    response_body = json.loads(request_handler.wfile.read().decode("utf-8"))
    return responses[0], response_body


def test_compare_endpoint_returns_result_schema() -> None:
    status_code, result = _post_compare(
        {
            "sequence": SAMPLE_PROTEIN,
            "profiles": ["balanced", "high_cai"],
            "scan_mode": "fast",
        }
    )

    assert status_code == 200
    assert list(result) == ["results"]
    assert [row["profile"] for row in result["results"]] == ["balanced", "high_cai"]
    for row in result["results"]:
        assert set(row) == {"profile", "cai", "gc_percent", "score", "sequence"}
        assert isinstance(row["cai"], float)
        assert isinstance(row["gc_percent"], float)
        assert isinstance(row["score"], float)
        assert row["sequence"].startswith("ATG")


def test_compare_endpoint_uses_default_profiles() -> None:
    h = _handler()

    status_code, result = h.handle_compare_request({"sequence": SAMPLE_PROTEIN})

    assert status_code == 200
    assert [row["profile"] for row in result["results"]] == [
        "balanced",
        "high_cai",
        "gc_target",
        "assembly_friendly",
    ]


def test_compare_endpoint_rejects_invalid_profile() -> None:
    h = _handler()

    status_code, result = h.handle_compare_request(
        {"sequence": SAMPLE_PROTEIN, "profiles": ["balanced", "bad_profile"]}
    )

    assert status_code == 400
    assert result == {"success": False, "error": "Invalid profile: bad_profile"}


def test_compare_endpoint_rejects_too_many_profiles() -> None:
    h = _handler()

    status_code, result = h.handle_compare_request(
        {
            "sequence": SAMPLE_PROTEIN,
            "profiles": [
                "balanced",
                "high_cai",
                "gc_target",
                "assembly_friendly",
                "balanced",
                "high_cai",
                "gc_target",
            ],
        }
    )

    assert status_code == 400
    assert result["error"] == "profiles must include at most 6 profiles"


def test_compare_endpoint_rejects_host_field() -> None:
    h = _handler()

    status_code, result = h.handle_compare_request(
        {"sequence": SAMPLE_PROTEIN, "profiles": ["balanced"], "host": "by2"}
    )

    assert status_code == 400
    assert result["error_code"] == "HOST_NOT_SUPPORTED_ON_ENDPOINT"


def test_compare_endpoint_rejects_host_profile_field() -> None:
    h = _handler()

    status_code, result = h.handle_compare_request(
        {"sequence": SAMPLE_PROTEIN, "profiles": ["balanced"], "host_profile": "by2"}
    )

    assert status_code == 400
    assert result["error_code"] == "HOST_NOT_SUPPORTED_ON_ENDPOINT"
