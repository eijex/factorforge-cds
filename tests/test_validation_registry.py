"""Tests for the runtime validation-check registry."""

from __future__ import annotations

import pytest

from factorforge.validation_registry import (
    LEGACY_VALIDATION_FIELD_MAP,
    PUBLIC_VALIDATION_BADGE_IDS,
    REVIEW_CONTRACT_CHECK_IDS,
    VALIDATION_CHECKS,
    VALIDATION_REGISTRY_VERSION,
    get_check,
    public_badge_checks,
)

_VALID_DOMAINS = {"integrity", "configured_constraint", "advisory_scan", "assembly_review"}
_REQUIRED_ENFORCEMENT_KEYS = {"bare_optimizer", "pipeline_default", "benchmark_scoring"}


def test_registry_version_is_set() -> None:
    assert VALIDATION_REGISTRY_VERSION == "1.0"


def test_registry_has_17_checks() -> None:
    assert len(VALIDATION_CHECKS) == 17


def test_review_contract_has_11_ids() -> None:
    assert len(REVIEW_CONTRACT_CHECK_IDS) == 11


def test_public_badge_ids_has_12_ids_including_gc() -> None:
    assert len(PUBLIC_VALIDATION_BADGE_IDS) == 12
    assert "global_gc_range" in PUBLIC_VALIDATION_BADGE_IDS
    assert set(REVIEW_CONTRACT_CHECK_IDS).issubset(set(PUBLIC_VALIDATION_BADGE_IDS))


def test_check_ids_are_unique() -> None:
    ids = [check["check_id"] for check in VALIDATION_CHECKS]
    assert len(ids) == len(set(ids))


def test_check_orders_are_unique() -> None:
    orders = [check["order"] for check in VALIDATION_CHECKS]
    assert len(orders) == len(set(orders))


def test_every_check_has_non_empty_display_name() -> None:
    for check in VALIDATION_CHECKS:
        assert isinstance(check["display_name"], str) and check["display_name"].strip()


def test_every_check_has_a_valid_domain() -> None:
    for check in VALIDATION_CHECKS:
        assert check["primary_domain"] in _VALID_DOMAINS


def test_every_check_has_all_three_enforcement_keys_and_no_others() -> None:
    for check in VALIDATION_CHECKS:
        assert set(check["enforcement_per_path"].keys()) == _REQUIRED_ENFORCEMENT_KEYS


def test_no_check_defines_a_web_ui_enforcement_path() -> None:
    for check in VALIDATION_CHECKS:
        assert "web_ui" not in check["enforcement_per_path"]
        assert "web" not in check["enforcement_per_path"]


def test_moclo_overhang_presentation_defaults_to_not_run() -> None:
    check = get_check("moclo_overhang")
    assert check["presentation"]["web"]["default_status_when_not_executed"] == "NOT_RUN"
    assert check["presentation"]["web"]["visible"] is True


def test_get_check_returns_the_matching_check() -> None:
    check = get_check("polya")
    assert check["check_id"] == "polya"
    assert check["primary_domain"] == "advisory_scan"


def test_get_check_raises_keyerror_for_unknown_id() -> None:
    with pytest.raises(KeyError):
        get_check("not_a_real_check")


def test_public_badge_checks_returns_12_sorted_by_order() -> None:
    checks = public_badge_checks()
    assert len(checks) == 12
    assert [check["check_id"] for check in checks] == sorted(
        (check["check_id"] for check in checks),
        key=lambda check_id: get_check(check_id)["order"],
    )
    orders = [check["order"] for check in checks]
    assert orders == sorted(orders)


def test_legacy_field_map_matches_job_142_fix_semantics() -> None:
    assert LEGACY_VALIDATION_FIELD_MAP == {
        "polya": "polya",
        "gc": "global_gc_range",
        "moclo": "restriction_sites",
    }
