"""Tests for the /api/optimize v1 contract helpers."""

from __future__ import annotations

from io import BytesIO
import json

import api.optimize as optimize_api
import pytest
from api.optimize import DEFAULT_GC_MAX, DEFAULT_GC_MIN, handler


def _handler() -> handler:
    return object.__new__(handler)


def _post_api(path: str, payload: dict) -> tuple[int, dict]:
    request_handler = _handler()
    body = json.dumps(payload).encode("utf-8")
    request_handler.path = path
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


def _post_optimize(payload: dict) -> tuple[int, dict]:
    return _post_api("/api/optimize", payload)


def test_public_api_rejects_codon_reference_override_fields() -> None:
    status_code, result = _post_optimize(
        {
            "sequence": "MSKGEELFTGVVPILVELD",
            "profile": "balanced",
            "codon_reference": "nbenthamiana_nbev11_hc_v2",
        }
    )

    assert status_code == 400
    assert result["error_code"] == "REFERENCE_OVERRIDE_NOT_SUPPORTED"
    assert result["forbidden_fields"] == ["codon_reference"]
    assert result["reference_policy"]["override_supported"] is False
    assert result["reference_policy"]["selectable_reference_ids"] == []
    assert (
        result["reference_policy"]["active_default_reference_id"]
        == "nbenthamiana_legacy_kazusa_sgn_v101"
    )


def test_compare_and_batch_reject_reference_override_fields() -> None:
    compare_status, compare_result = _post_api(
        "/api/optimize/compare",
        {
            "sequence": "MSKGEELFTGVVPILVELD",
            "profiles": ["balanced"],
            "codon_table_id": "nbenthamiana_nbev11_hc_v2",
        },
    )
    batch_status, batch_result = _post_api(
        "/api/optimize/batch",
        {
            "sequences": [{"id": "seq1", "sequence": "MSKGEELFTGVVPILVELD"}],
            "profile": "balanced",
            "codon_table_path": "src/factorforge/data/profiles/nbev11_cds_hc_derived_codons.json",
        },
    )

    assert compare_status == 400
    assert compare_result["error_code"] == "REFERENCE_OVERRIDE_NOT_SUPPORTED"
    assert compare_result["forbidden_fields"] == ["codon_table_id"]
    assert batch_status == 400
    assert batch_result["error_code"] == "REFERENCE_OVERRIDE_NOT_SUPPORTED"
    assert batch_result["forbidden_fields"] == ["codon_table_path"]


def test_public_reference_policy_metadata_does_not_advertise_disabled_selectors() -> None:
    metadata = optimize_api._public_reference_policy_metadata()

    assert metadata["override_supported"] is False
    assert metadata["selectable_reference_ids"] == []
    assert metadata["active_default"] == {
        "reference_id": "nbenthamiana_legacy_kazusa_sgn_v101",
        "tier": "production_enabled",
        "activation_status": "enabled",
        "claim_boundary": metadata["active_default"]["claim_boundary"],
    }
    assert metadata["activation_status_counts"]["enabled"] == 1
    assert metadata["activation_status_counts"]["disabled"] >= 1
    assert metadata["activation_status_counts"]["research_only"] >= 1
    assert metadata["disabled_or_research_only_count"] >= 3
    assert "codon_reference" in metadata["forbidden_override_fields"]


def test_get_optimize_exposes_public_reference_policy_metadata() -> None:
    data = _get_optimize()
    policy = data["reference_policy"]

    assert policy["override_supported"] is False
    assert policy["selectable_reference_ids"] == []
    assert policy["active_default_reference_id"] == "nbenthamiana_legacy_kazusa_sgn_v101"
    assert policy["active_default"]["activation_status"] == "enabled"


def test_parse_constraints_defaults() -> None:
    h = _handler()

    # No host ⇒ defaults to nbenthamiana, whose default band is the
    # _analysis/025 composition anchor (Job 168 / v3.3.0), not the
    # module-level DEFAULT_GC_MIN/MAX (which is now only the fallback for
    # hosts without their own analysis, e.g. by2/ntabacum).
    assert h.parse_constraints({}) == optimize_api._default_gc_constraints("nbenthamiana")


def test_parse_constraints_defaults_for_non_nbenthamiana_host_unchanged() -> None:
    h = _handler()

    # by2/ntabacum has no host-specific analysis (Job 168 scope: nbenthamiana
    # only) — it must keep the pre-v3.3.0 global default unchanged.
    assert h.parse_constraints({}, host="ntabacum") == {
        "gc_min": DEFAULT_GC_MIN,
        "gc_max": DEFAULT_GC_MAX,
    }


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
    assert result["engine_versions"]["product"] == "3.2.8"
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


def test_balanced_api_response_exposes_gc_target_observation() -> None:
    h = _handler()

    result = h.optimize_sequence(
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG",
        "balanced",
        False,
        False,
        False,
        objective=None,
        return_candidates=False,
        constraints={"gc_min": 40.0, "gc_max": 47.0},
    )

    metrics = result["metrics"]
    assert metrics["requested_gc_min_percent"] == 40.0
    assert metrics["requested_gc_max_percent"] == 47.0
    assert metrics["gc_target_reached"] == (40.0 <= metrics["gc_percent"] <= 47.0)


def test_non_balanced_api_response_omits_gc_target_observation() -> None:
    h = _handler()

    result = h.optimize_sequence(
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG",
        "gc_target",
        False,
        False,
        False,
        objective=None,
        return_candidates=False,
        constraints={"gc_min": 40.0, "gc_max": 47.0},
    )

    metrics = result["metrics"]
    assert "gc_target_reached" not in metrics
    assert "requested_gc_min_percent" not in metrics
    assert "requested_gc_max_percent" not in metrics


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
    optimized_sequence = "ATGCTACTACTA"

    # Whether a run of CTA codons counts as "rare" depends on the active
    # production codon-usage table's relative-adaptiveness weight for CTA
    # (Job 168 / v3.3.0 migrated nbenthamiana's default — see
    # data/reference/active_codon_reference.json). Derive the expected count
    # from the live rule engine instead of hardcoding it, so this test keeps
    # verifying that add_design_package_fields actually wires through to
    # RuleEngine.scan_rare_codon_runs rather than re-encoding a stale
    # rarity-classification assumption that drifts with the active table.
    from factorforge.engines.profile.rules.rule_engine import RuleEngine

    expected_clusters = len(RuleEngine().scan_rare_codon_runs(optimized_sequence))

    result = h.add_design_package_fields(
        response={
            "success": True,
            "optimized_sequence": optimized_sequence,
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

    assert result["constraint_report"]["codon_rarity_clusters"] == expected_clusters


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


# --- BY-2 strategy/host compatibility guard + disclosure ---


def test_explicit_feasibility_best_rejected_for_by2() -> None:
    status_code, result = _post_optimize(
        {"sequence": "MSKGEELFTGVVPILVELD", "objective": "feasibility_best", "host": "by2"}
    )

    assert status_code == 400
    assert result["success"] is False
    assert isinstance(result["error"], str)
    assert "feasibility_best" in result["error"]
    assert result["error_code"] == "UNSUPPORTED_STRATEGY_HOST_COMBINATION"
    assert result["requested_host"] == "ntabacum"
    assert result["requested_strategy"] == "feasibility_best"


def test_explicit_high_cai_rejected_for_by2() -> None:
    status_code, result = _post_optimize(
        {"sequence": "MSKGEELFTGVVPILVELD", "profile": "high_cai", "host": "by2"}
    )

    assert status_code == 400
    assert result["success"] is False
    assert isinstance(result["error"], str)
    assert "high_cai" in result["error"]
    assert result["error_code"] == "UNSUPPORTED_STRATEGY_HOST_COMBINATION"
    assert result["requested_host"] == "ntabacum"
    assert result["requested_strategy"] == "high_cai"


def test_explicit_feasibility_best_unaffected_for_nbenthamiana() -> None:
    status_code, result = _post_optimize(
        {"sequence": "MSKGEELFTGVVPILVELD", "objective": "feasibility_best", "host": "nbenthamiana"}
    )

    assert status_code == 200
    assert result["success"] is True
    assert "error_code" not in result


def test_explicit_high_cai_unaffected_for_nbenthamiana() -> None:
    status_code, result = _post_optimize(
        {"sequence": "MSKGEELFTGVVPILVELD", "profile": "high_cai", "host": "nbenthamiana"}
    )

    assert status_code == 200
    assert result["success"] is True
    assert "error_code" not in result


def test_implicit_by2_resolves_to_balanced_with_disclosure() -> None:
    status_code, result = _post_optimize({"sequence": "MSKGEELFTGVVPILVELD", "host": "by2"})

    assert status_code == 200
    assert result["success"] is True
    assert result["profile"] == "balanced"
    assert result["requested_strategy"] == "feasibility_best"
    assert result["resolved_strategy"] == "balanced"
    assert "not available for this host" in result["resolution_reason"]


def test_explicit_balanced_for_by2_is_host_aware_and_unrejected() -> None:
    status_code, result = _post_optimize(
        {"sequence": "MSKGEELFTGVVPILVELD", "profile": "balanced", "host": "by2"}
    )

    assert status_code == 200
    assert result["success"] is True
    assert "error_code" not in result
    assert "resolved_strategy" not in result


def test_invalid_host_error_contract_unchanged() -> None:
    status_code, result = _post_optimize(
        {"sequence": "MSKGEELFTGVVPILVELD", "profile": "balanced", "host": "not_a_host"}
    )

    assert status_code == 400
    assert result["success"] is False
    assert isinstance(result["error"], str)
    assert "Invalid host" in result["error"]


# ── host_metadata.gc_range web sync (Job 168 / v3.3.0, STEP 4.1) ─────────────
# web/js/app.js must never hardcode a GC reference band — it reads
# host_metadata[host].gc_range from this GET response at runtime. These tests
# guard the response itself, plus the offline-only fallback literal in app.js,
# against drifting away from resolve_host_gc_range (the single source of truth).


def _get_optimize() -> dict:
    request_handler = _handler()
    request_handler.headers = {}
    request_handler.wfile = BytesIO()
    responses: list[int] = []
    headers: list[tuple[str, str]] = []
    request_handler.send_response = responses.append
    request_handler.send_header = lambda key, value: headers.append((key, value))
    request_handler.end_headers = lambda: None

    handler.do_GET(request_handler)

    request_handler.wfile.seek(0)
    return json.loads(request_handler.wfile.read().decode("utf-8"))


def test_host_metadata_gc_range_matches_resolve_host_gc_range() -> None:
    from factorforge.engines.profile.scoring import resolve_host_gc_range

    data = _get_optimize()
    nb_min, nb_max = resolve_host_gc_range("nbenthamiana")
    by2_min, by2_max = resolve_host_gc_range("ntabacum")

    assert data["host_metadata"]["nbenthamiana"]["gc_range"] == {
        "gc_min": nb_min,
        "gc_max": nb_max,
    }
    assert data["host_metadata"]["by2"]["gc_range"] == {"gc_min": by2_min, "gc_max": by2_max}


def test_appjs_offline_gc_fallback_matches_resolve_host_gc_range() -> None:
    """app.js's OFFLINE_GC_RANGES is a same-values offline/dev fallback only
    (used when GET /api/optimize is unreachable) — it must mirror the live
    resolver, not drift into its own competing source of truth."""
    import re
    from pathlib import Path

    from factorforge.engines.profile.scoring import resolve_host_gc_range

    app_js = (Path(__file__).resolve().parents[2] / "web" / "js" / "app.js").read_text(
        encoding="utf-8"
    )
    match = re.search(
        r"const OFFLINE_GC_RANGES = \{(.*?)\n\};", app_js, re.DOTALL
    )
    assert match is not None, "OFFLINE_GC_RANGES literal not found in web/js/app.js"
    block = match.group(1)

    for host_id, internal_host in (("nbenthamiana", "nbenthamiana"), ("by2", "ntabacum")):
        host_match = re.search(
            rf"{host_id}:\s*\{{\s*gc_min:\s*([\d.]+),\s*gc_max:\s*([\d.]+)",
            block,
        )
        assert host_match is not None, f"OFFLINE_GC_RANGES.{host_id} not found"
        js_min, js_max = float(host_match.group(1)), float(host_match.group(2))
        py_min, py_max = resolve_host_gc_range(internal_host)
        assert (js_min, js_max) == (py_min, py_max), (
            f"OFFLINE_GC_RANGES.{host_id} ({js_min}, {js_max}) != "
            f"resolve_host_gc_range({internal_host!r}) ({py_min}, {py_max})"
        )
