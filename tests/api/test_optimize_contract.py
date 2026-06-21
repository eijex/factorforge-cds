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
    assert result["engine_versions"]["product"] == "3.2.3"
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
    assert result["recommended_candidate"]["cai"] == result["metrics"]["cai"]
    assert result["recommended_candidate"]["cai_reference"] == "profile_golden_set"
    assert result["metrics"]["cai_reference"] == "profile_golden_set"
    assert result["metrics"]["general_cai"] == result["recommended_candidate"]["general_cai"]
    assert isinstance(result["recommended_candidate"]["general_cai"], float)
    assert isinstance(result["recommended_candidate"]["assembly_pass"], bool)
    assert isinstance(result["recommended_candidate"]["forbidden_type_iis_site_count"], int)
    expected_moclo = "PASS" if result["recommended_candidate"]["assembly_pass"] else "WARNING"
    assert result["validation"]["moclo"] == expected_moclo


def test_candidate_with_type_iis_conflict_is_not_reported_as_pass() -> None:
    h = _handler()
    table = optimize_api.load_codon_usage_table()

    candidate = h.build_candidate(
        candidate_id="test",
        label="Test",
        dna_sequence="ATGGGTCTCGCT",
        codon_weights=table.codon_weights,
        recommendation_reason="test fixture",
    )

    assert candidate["assembly_pass"] is False
    assert candidate["forbidden_type_iis_site_count"] == 1
    assert candidate["validator_status"] == "warning"


def test_build_candidate_includes_canonical_checks() -> None:
    h = _handler()
    table = optimize_api.load_codon_usage_table()

    candidate = h.build_candidate(
        candidate_id="test",
        label="Test",
        dna_sequence="ATGGCTGGCGGCGGCAGCGGCAGCAATAAAGGCGGCAGCGGCAGCGGCGGCAGCGGCGGC",
        codon_weights=table.codon_weights,
        recommendation_reason="test fixture",
        constraints={"gc_min": 40.0, "gc_max": 70.0},
    )

    assert "checks" in candidate
    assert candidate["checks"]["polya"]["status"] == "WARNING"
    assert candidate["checks"]["polya"]["finding_count"] >= 1
    assert candidate["checks"]["moclo_overhang"]["status"] == "NOT_RUN"


def test_feasibility_best_candidates_each_carry_checks() -> None:
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

    assert "checks" in result["recommended_candidate"]
    for candidate in result["candidates"]:
        assert "checks" in candidate


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


def test_get_health_check_exposes_validation_registry() -> None:
    h = _handler()
    responses: list[int] = []
    headers: list[tuple[str, str]] = []
    sent: dict[str, object] = {}
    h.path = "/api/optimize"
    h.send_response = responses.append
    h.send_header = lambda key, value: headers.append((key, value))
    h.end_headers = lambda: None
    h.wfile = BytesIO()

    handler.do_GET(h)

    h.wfile.seek(0)
    body = json.loads(h.wfile.read().decode("utf-8"))
    assert body["validation_registry_version"] == "1.0"
    assert body["validation_report_schema_version"] == "1.0"
    assert len(body["validation_checks"]) == 12
    assert body["validation_checks"][0]["check_id"] == "global_gc_range"


def test_legacy_profile_response_gains_additive_validation_checks() -> None:
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

    # Existing legacy fields untouched (AC9):
    assert result["validation"]["polya"] in {"PASS", "WARNING"}
    assert result["validation"]["moclo"] in {"PASS", "WARNING"}
    assert result["validation"]["gc"] in {"PASS", "WARNING"}
    # New additive fields:
    assert result["validation"]["schema_version"] == "1.0"
    assert "checks" in result["validation"]
    assert result["metadata"]["validation_registry_version"] == "1.0"


def test_feasibility_best_response_gains_additive_validation_checks() -> None:
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

    # Existing fields untouched (AC9):
    assert result["validation"] == {
        "input_type": "protein",
        "sequence_length": 35,
        "host_profile": "nbenthamiana",
    }
    # New additive top-level fields, alongside the unchanged dict above —
    # added as siblings, not merged into the existing "validation" dict,
    # since that dict's exact shape is pinned by AC9.
    assert result["validation_report"]["schema_version"] == "1.0"
    assert "checks" in result["validation_report"]
    assert result["metadata"]["validation_registry_version"] == "1.0"


def test_custom_restriction_site_rebuild_uses_request_host_and_constraints() -> None:
    """Rebuilt candidate checks must use the request's host/constraints, not defaults.

    GAGGAG is present in the gc_target-optimized DNA for this protein under the
    ntabacum host, so domestication actually fires and update_candidate_sequence
    rebuilds the candidate via build_candidate. The GC band 70-80% is chosen so
    that the rebuilt candidate's ~59% GC is a WARNING under the *custom* band but
    would be a PASS under the default 55-65% band — this discriminates between
    "used the real constraints" and "silently fell back to defaults".
    """
    h = _handler()

    result = h.optimize_sequence(
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG",
        "gc_target",
        False,
        False,
        False,
        objective=None,
        host="ntabacum",
        return_candidates=True,
        constraints={"gc_min": 70.0, "gc_max": 80.0},
        custom_restriction_sites=[{"name": "TestSite", "sequence": "GAGGAG"}],
    )

    # Domestication must have actually fired for this assertion to be meaningful.
    assert result["custom_restriction_sites"]["removed"]

    recommended = result["recommended_candidate"]
    # The rebuilt candidate's GC% sits inside the *default* 55-65% band (would be
    # PASS if the rebuild silently fell back to defaults) but outside the custom
    # 70-80% band actually requested — proving the real constraints were used.
    assert DEFAULT_GC_MIN <= recommended["gc_percent"] <= DEFAULT_GC_MAX
    assert recommended["checks"]["global_gc_range"]["status"] == "WARNING"
    for candidate in result["candidates"]:
        assert candidate["checks"]["global_gc_range"]["status"] == "WARNING"
