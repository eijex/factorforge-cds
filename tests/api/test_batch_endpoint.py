"""Tests for the /api/optimize/batch endpoint."""

from __future__ import annotations

from io import BytesIO
import json

from api.optimize import handler


SAMPLE_PROTEIN = "MSKGEELFTGVVPILVELD"


def _handler() -> handler:
    return object.__new__(handler)


def _post_batch(payload: dict) -> tuple[int, dict]:
    request_handler = _handler()
    body = json.dumps(payload).encode("utf-8")
    request_handler.path = "/api/optimize/batch"
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


def test_batch_endpoint_returns_result_schema() -> None:
    status_code, result = _post_batch(
        {
            "sequences": [
                {"id": "GFP", "sequence": SAMPLE_PROTEIN},
                {"id": "mCherry", "sequence": "MVSSKGEELFTGVVPILVELD"},
            ],
            "profile": "balanced",
            "scan_mode": "fast",
        }
    )

    assert status_code == 200
    assert result["count"] == 2
    assert result["profile"] == "balanced"
    assert [row["id"] for row in result["results"]] == ["GFP", "mCherry"]
    for row in result["results"]:
        assert set(row) == {"id", "sequence", "cai", "gc_percent", "score", "violations"}
        assert row["sequence"].startswith("ATG")
        assert isinstance(row["cai"], float)
        assert isinstance(row["gc_percent"], float)
        assert isinstance(row["score"], float)
        assert isinstance(row["violations"], int)


def test_batch_endpoint_generates_missing_ids() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request(
        {"sequences": [{"sequence": SAMPLE_PROTEIN}, {"sequence": SAMPLE_PROTEIN}]}
    )

    assert status_code == 200
    assert [row["id"] for row in result["results"]] == ["seq_1", "seq_2"]
    assert result["profile"] == "balanced"


def test_batch_endpoint_rejects_21_sequences() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request(
        {"sequences": [{"sequence": SAMPLE_PROTEIN} for _ in range(21)]}
    )

    assert status_code == 400
    assert result == {"success": False, "error": "Batch limit is 20 sequences"}


def test_batch_endpoint_rejects_empty_sequences() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request({"sequences": []})

    assert status_code == 400
    assert result == {
        "success": False,
        "error": "sequences is required and must be non-empty",
    }


def test_batch_endpoint_rejects_missing_sequences() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request({})

    assert status_code == 400
    assert result == {
        "success": False,
        "error": "sequences is required and must be non-empty",
    }


def test_batch_endpoint_rejects_invalid_profile() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request(
        {"sequences": [{"sequence": SAMPLE_PROTEIN}], "profile": "bad_profile"}
    )

    assert status_code == 400
    assert result["error"].startswith("Invalid profile. Must be one of:")


def test_batch_endpoint_rejects_host_field() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request(
        {"sequences": [{"sequence": SAMPLE_PROTEIN}], "host": "by2"}
    )

    assert status_code == 400
    assert result["error_code"] == "HOST_NOT_SUPPORTED_ON_ENDPOINT"


def test_batch_endpoint_rejects_host_profile_field() -> None:
    h = _handler()

    status_code, result = h.handle_batch_request(
        {"sequences": [{"sequence": SAMPLE_PROTEIN}], "host_profile": "by2"}
    )

    assert status_code == 400
    assert result["error_code"] == "HOST_NOT_SUPPORTED_ON_ENDPOINT"
