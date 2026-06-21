"""Tests for the canonical, check_id-keyed validation report builder."""

from __future__ import annotations

from factorforge.validation_report import (
    VALIDATION_REPORT_SCHEMA_VERSION,
    build_validation_report,
)

_DEFAULT_CONSTRAINTS = {"gc_min": 55.0, "gc_max": 65.0}

# 60bp CDS, in-frame, GC% in [55,65], with a deliberate AATAAA PolyA motif and
# no Type IIS / homopolymer / overhang side effects — isolates the PolyA signal.
_POLYA_POSITIVE_CDS = (
    "ATGGCTGGCGGCGGCAGCGGCAGCAATAAAGGCGGCAGCGGCAGCGGCGGCAGCGGCGGC"
)
assert len(_POLYA_POSITIVE_CDS) % 3 == 0

_CLEAN_CDS = "ATGGCTGGCGGCGGCAGCGGCAGCGGCGGCAGCGGCAGCGGCGGCAGCGGCGGCAGCGGC"
assert len(_CLEAN_CDS) % 3 == 0

_TYPE_IIS_POSITIVE_CDS = "ATGGGTCTCGCTGGCGGCAGCGGCGGCAGCGGCAGCGGCGGCAGCGGCAGC"


def test_schema_version_is_set() -> None:
    assert VALIDATION_REPORT_SCHEMA_VERSION == "1.0"


def test_report_contains_all_nine_advisory_scan_keys() -> None:
    report = build_validation_report(
        _CLEAN_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS
    )
    for check_id in (
        "polya", "are", "at_runs", "homopolymers", "repeats",
        "gc_extremes", "splice_sites", "dinucleotides", "rare_codon_runs",
    ):
        assert check_id in report["checks"]


def test_report_contains_gc_and_restriction_sites_keys() -> None:
    report = build_validation_report(
        _CLEAN_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS
    )
    assert "global_gc_range" in report["checks"]
    assert "restriction_sites" in report["checks"]


def test_polya_positive_fixture_produces_warning_with_finding(value_drift_regression=None) -> None:
    """AC4 / AC12: scanners must run against the exact final_cds passed in, and a
    real PolyA motif in that sequence must surface as a non-zero, non-PASS result."""
    report = build_validation_report(
        _POLYA_POSITIVE_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS
    )
    polya = report["checks"]["polya"]
    assert polya["status"] == "WARNING"
    assert polya["finding_count"] >= 1
    assert polya["executed"] is True


def test_clean_sequence_polya_passes_with_zero_findings() -> None:
    report = build_validation_report(
        _CLEAN_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS
    )
    polya = report["checks"]["polya"]
    assert polya["status"] == "PASS"
    assert polya["finding_count"] == 0


def test_global_gc_range_pass_and_warning() -> None:
    passing = build_validation_report(_CLEAN_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS)
    assert passing["checks"]["global_gc_range"]["status"] == "PASS"

    failing = build_validation_report(_CLEAN_CDS, gc_percent=10.0, constraints=_DEFAULT_CONSTRAINTS)
    assert failing["checks"]["global_gc_range"]["status"] == "WARNING"


def test_restriction_sites_detects_type_iis_motif() -> None:
    report = build_validation_report(
        _TYPE_IIS_POSITIVE_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS
    )
    restriction = report["checks"]["restriction_sites"]
    assert restriction["status"] == "WARNING"
    assert restriction["finding_count"] >= 1


def test_moclo_overhang_defaults_to_not_run_when_not_requested() -> None:
    report = build_validation_report(
        _CLEAN_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS, moclo_requested=False
    )
    moclo = report["checks"]["moclo_overhang"]
    assert moclo["status"] == "NOT_RUN"
    assert moclo["finding_count"] is None
    assert moclo["executed"] is False
    assert moclo["reason"] == "construct_template_not_requested"


def test_moclo_overhang_runs_and_flags_internal_collision_when_requested() -> None:
    cds_with_internal_overhang = "ATGAATGGGCGGCAGCGGCAGCGGCGGCAGCGGCAGCGGCGGCAGCGGC"
    report = build_validation_report(
        cds_with_internal_overhang,
        gc_percent=60.0,
        constraints=_DEFAULT_CONSTRAINTS,
        moclo_requested=True,
    )
    moclo = report["checks"]["moclo_overhang"]
    assert moclo["executed"] is True
    assert moclo["status"] == "WARNING"
    assert moclo["finding_count"] >= 1


def test_type_iis_and_moclo_are_independent_checks() -> None:
    """AC7: Type IIS and MoClo retain different enforcement semantics — a
    Type-IIS-positive sequence must not implicitly mark MoClo as executed."""
    report = build_validation_report(
        _TYPE_IIS_POSITIVE_CDS,
        gc_percent=60.0,
        constraints=_DEFAULT_CONSTRAINTS,
        moclo_requested=False,
    )
    assert report["checks"]["restriction_sites"]["status"] == "WARNING"
    assert report["checks"]["moclo_overhang"]["status"] == "NOT_RUN"


def test_all_check_statuses_are_in_the_canonical_enum() -> None:
    valid_statuses = {"PASS", "WARNING", "FAIL", "NOT_RUN", "NOT_APPLICABLE", "ERROR"}
    report = build_validation_report(
        _POLYA_POSITIVE_CDS, gc_percent=60.0, constraints=_DEFAULT_CONSTRAINTS, moclo_requested=True
    )
    for check_id, result in report["checks"].items():
        assert result["status"] in valid_statuses, f"{check_id} has invalid status {result['status']}"
