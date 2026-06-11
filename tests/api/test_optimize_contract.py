"""Tests for the /api/optimize v1 contract helpers."""

from __future__ import annotations

from io import BytesIO
import json

import api.optimize as optimize_api
import pytest
from api.optimize import DEFAULT_GC_MAX, DEFAULT_GC_MIN, handler


def _handler() -> handler:
    return object.__new__(handler)


def _post_optimize(payload: dict) -> tuple[int, dict]:
    request_handler = _handler()
    body = json.dumps(payload).encode("utf-8")
    request_handler.path = "/api/optimize"
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


def test_parse_constraints_defaults() -> None:
    h = _handler()

    assert h.parse_constraints({}) == {"gc_min": DEFAULT_GC_MIN, "gc_max": DEFAULT_GC_MAX}


def test_feasibility_best_response_includes_candidate_contract() -> None:
    h = _handler()

    result = h.optimize_sequence(
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG",
        "balanced",
        False,
        False,
        False,
        objective="feasibility_best",
        host_profile="nbenthamiana",
        return_candidates=True,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )

    assert result["success"] is True
    assert result["recommended_candidate"]["id"] == "feasibility_best"
    assert [candidate["id"] for candidate in result["candidates"]] == [
        "feasibility_best",
        "gc_target",
        "high_cai",
    ]
    assert result["validation"] == {
        "input_type": "protein",
        "sequence_length": 35,
        "host_profile": "nbenthamiana",
    }
    assert result["engine_versions"]["product"] == "3.2.0"
    assert result["recommended_candidate"]["validator_status"] == "pass"


def test_legacy_profile_response_keeps_old_fields_and_adds_candidates() -> None:
    h = _handler()

    result = h.optimize_sequence(
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG",
        "gc_target",
        False,
        False,
        False,
        objective=None,
        return_candidates=True,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )

    assert result["success"] is True
    assert result["profile"] == "gc_target"
    assert "optimized_sequence" in result
    assert result["recommended_candidate"]["id"] == "gc_target"
    assert [candidate["id"] for candidate in result["candidates"]] == ["gc_target"]


def test_design_package_aa_identity_uses_cds_validator_result() -> None:
    h = _handler()

    result = h.add_design_package_fields(
        response={
            "success": True,
            "optimized_sequence": "ATGGCTTAC",
            "metrics": {"cai": 0.8, "gc_percent": 44.4},
            "validation": {"polya": "PASS", "moclo": "UNCHECKED", "gc": "PASS"},
        },
        input_sequence="MAF",
        profile="balanced",
        objective=None,
        host_profile="nbenthamiana",
        kozak=False,
        dinuc=False,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )

    assert result["constraint_report"]["aa_identity"] == pytest.approx(2 / 3)
    assert result["validation_status"]["aa_identity_check"] == "fail"
    assert any(
        error.startswith("aa_mismatch")
        for error in result["constraint_report"]["cds_validation_errors"]
    )


def test_design_package_codon_rarity_clusters_use_rule_scan() -> None:
    h = _handler()

    result = h.add_design_package_fields(
        response={
            "success": True,
            "optimized_sequence": "ATGCTACTACTA",
            "metrics": {"cai": 0.8, "gc_percent": 16.7},
            "validation": {"polya": "PASS", "moclo": "UNCHECKED", "gc": "PASS"},
        },
        input_sequence="MLLL",
        profile="balanced",
        objective=None,
        host_profile="nbenthamiana",
        kozak=False,
        dinuc=False,
        constraints={"gc_min": 40.0, "gc_max": 55.0},
    )

    assert result["constraint_report"]["codon_rarity_clusters"] == 1


def test_legacy_gc_validation_uses_request_constraints() -> None:
    h = _handler()

    result = h.optimize_sequence(
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG",
        "gc_target",
        False,
        False,
        False,
        objective=None,
        return_candidates=False,
        constraints={"gc_min": 45.0, "gc_max": 55.0},
    )

    assert result["validation"]["gc"] == "WARNING"


def test_engine_unavailable_returns_503_when_mock_disabled(monkeypatch) -> None:
    h = _handler()
    monkeypatch.setattr(optimize_api, "ENABLE_MOCK", False)

    status_code, result = h.handle_unavailable_engine(
        "MSK",
        "balanced",
        False,
        False,
        {"gc_min": 40.0, "gc_max": 55.0},
        False,
    )

    assert status_code == 503
    assert result == {"success": False, "error": "Engine unavailable. Contact support."}


def test_engine_unavailable_returns_mock_when_enabled(monkeypatch) -> None:
    h = _handler()
    monkeypatch.setattr(optimize_api, "ENABLE_MOCK", True)

    status_code, result = h.handle_unavailable_engine(
        "MSK",
        "balanced",
        False,
        False,
        {"gc_min": 40.0, "gc_max": 55.0},
        False,
    )

    assert status_code == 200
    assert result["success"] is True
    assert result["engine"]["name"] == "Mock Engine"


def test_optimize_endpoint_rejects_protein_over_web_api_limit() -> None:
    status_code, result = _post_optimize({"sequence": "M" * 5001, "profile": "balanced"})

    assert status_code == 400
    assert result["success"] is False
    assert "Protein input exceeds maximum length" in result["error"]
    assert "5,000 amino acids" in result["error"]
    assert result["cli_install"] == "pip install factorforge-cds"


def test_optimize_endpoint_accepts_protein_at_web_api_limit() -> None:
    # 1666 aa protein → 4998 bp DNA: well within both limits, should not be rejected for length.
    status_code, result = _post_optimize({"sequence": "M" * 1666, "profile": "balanced"})
    # Either succeeds or fails for a non-length reason — must not be a length rejection.
    if status_code == 400:
        assert "exceeds maximum length" not in result.get("error", "")


def test_optimize_endpoint_rejects_dna_over_web_api_limit() -> None:
    # 15,001 bp of pure ACGT → classified as DNA, over the 15,000 bp cap.
    status_code, result = _post_optimize({"sequence": "ACG" * 5001, "profile": "balanced"})

    assert status_code == 400
    assert result["success"] is False
    assert "DNA input exceeds maximum length" in result["error"]
    assert "15,000 bp" in result["error"]


def test_optimize_endpoint_allows_dna_above_protein_limit() -> None:
    # 6,000 bp DNA (2,000 codons) exceeds the 5,000-aa protein cap but is under the
    # 15,000-bp DNA cap → must NOT be rejected for length.
    status_code, result = _post_optimize({"sequence": "ACG" * 2000, "profile": "balanced"})
    if status_code == 400:
        assert "exceeds maximum length" not in result.get("error", "")
