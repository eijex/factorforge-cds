"""Canonical, check_id-keyed validation report builder.

Converts raw RuleEngine / Domesticator / ConstructBuilder outputs into the
stable result envelope described in
docs/superpowers/specs/2026-06-21-factorforge-validation-registry-design.md
(eijex-workspace, private) §4. Every advisory scanner and the restriction-site
scan run fresh against `final_cds` — the exact sequence the caller is about to
return — never a stale or intermediate candidate (spec §4.3).
"""

from __future__ import annotations

from typing import Any

from factorforge.engines.profile.construct_builder import ConstructBuilder
from factorforge.engines.profile.rules.domesticator import Domesticator
from factorforge.engines.profile.rules.rule_engine import RuleEngine

VALIDATION_REPORT_SCHEMA_VERSION = "1.0"


def _scan_status(findings: list[Any]) -> dict[str, Any]:
    return {
        "status": "PASS" if not findings else "WARNING",
        "finding_count": len(findings),
        "executed": True,
        "reason": None,
    }


def build_validation_report(
    final_cds: str,
    *,
    gc_percent: float,
    constraints: dict[str, float],
    rule_engine: RuleEngine | None = None,
    domesticator: Domesticator | None = None,
    construct_builder: ConstructBuilder | None = None,
    moclo_requested: bool = False,
) -> dict[str, Any]:
    """Build the canonical validation report for `final_cds`.

    Args:
        final_cds: The exact DNA sequence the caller will return to the user.
        gc_percent: Global GC% of `final_cds` (caller already computes this).
        constraints: `{"gc_min": float, "gc_max": float}` for the request.
        rule_engine: Optional pre-configured RuleEngine (e.g. host-specific).
            Defaults to a host-default RuleEngine().
        domesticator: Optional pre-configured Domesticator. Defaults to Domesticator().
        construct_builder: Optional pre-configured ConstructBuilder. Defaults to
            ConstructBuilder() (no template_dir needed for overhang collision checks).
        moclo_requested: Whether construct/MoClo validation was explicitly
            requested for this request (e.g. use_template=True). When False,
            moclo_overhang is reported as NOT_RUN, never PASS.

    Returns:
        {"schema_version": "1.0", "checks": {check_id: {status, finding_count,
        executed, reason}, ...}} covering all 9 advisory scanners plus
        global_gc_range, restriction_sites, and moclo_overhang (12 keys).
    """
    rule_engine = rule_engine or RuleEngine()
    domesticator = domesticator or Domesticator()

    scan_results = rule_engine.scan_all(final_cds, mode="full")
    checks: dict[str, Any] = {
        scan_id: _scan_status(findings) for scan_id, findings in scan_results.items()
    }

    gc_in_range = constraints["gc_min"] <= gc_percent <= constraints["gc_max"]
    checks["global_gc_range"] = {
        "status": "PASS" if gc_in_range else "WARNING",
        "finding_count": 0 if gc_in_range else 1,
        "executed": True,
        "reason": None,
    }

    type_iis_sites = domesticator.scan_restriction_sites(final_cds, "golden_gate")
    checks["restriction_sites"] = _scan_status(type_iis_sites)

    if moclo_requested:
        builder = construct_builder or ConstructBuilder()
        collisions = builder.check_internal_overhang_collisions(final_cds)
        checks["moclo_overhang"] = _scan_status(collisions)
    else:
        checks["moclo_overhang"] = {
            "status": "NOT_RUN",
            "finding_count": None,
            "executed": False,
            "reason": "construct_template_not_requested",
        }

    return {"schema_version": VALIDATION_REPORT_SCHEMA_VERSION, "checks": checks}
