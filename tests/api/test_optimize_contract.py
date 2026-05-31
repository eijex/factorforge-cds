"""Tests for the /api/optimize v1 contract helpers."""

from __future__ import annotations

import api.optimize as optimize_api
from api.optimize import DEFAULT_GC_MAX, DEFAULT_GC_MIN, handler


def _handler() -> handler:
    return object.__new__(handler)


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
    assert result["engine_versions"]["product"] == "3.1.7"
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
