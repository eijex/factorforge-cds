"""Runtime validation-check registry — single source of truth for check_id,
display metadata, and per-execution-path enforcement.

Frozen reference: validation_contract_v1.yaml, factorforge commit
4a8be9f053797d5f54154afcbda732eaaf79f8ae (v3.2.3). This module is the
runtime source of truth; the contract file is a pinned manuscript snapshot
of what this registry's predecessor behavior was.
"""

from __future__ import annotations

from typing import Any

VALIDATION_REGISTRY_VERSION = "1.0"

VALIDATION_CHECKS: tuple[dict[str, Any], ...] = (
    # --- Domain B: Configured Constraints -----------------------------------
    {
        "check_id": "global_gc_range",
        "display_name": "Global GC% target band",
        "primary_domain": "configured_constraint",
        "order": 1,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "metric_only",
            "pipeline_default": "metric_only",
            "benchmark_scoring": "component_of_multi_constraint_pass_soft",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": (
            "Configured optimization target, not a pass/fail biological-risk scanner. "
            "Only gates in the benchmark scoring path."
        ),
    },
    # --- Domain C: Advisory Sequence Scans (RuleEngine, 9 scanners) ---------
    {
        "check_id": "polya",
        "display_name": "PolyA-like motifs",
        "primary_domain": "advisory_scan",
        "order": 2,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Advisory finding, not a guarantee of expression impact.",
    },
    {
        "check_id": "are",
        "display_name": "AU-rich elements (ARE)",
        "primary_domain": "advisory_scan",
        "order": 3,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Advisory finding.",
    },
    {
        "check_id": "at_runs",
        "display_name": "AT-rich runs",
        "primary_domain": "advisory_scan",
        "order": 4,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Advisory finding.",
    },
    {
        "check_id": "homopolymers",
        "display_name": "Homopolymer runs (synthesis risk)",
        "primary_domain": "advisory_scan",
        "order": 5,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": (
            "Synthesis-difficulty advisory, distinct from the separate expression-stability "
            "homopolymer threshold used in the unrelated/archived utils/validation.py path."
        ),
    },
    {
        "check_id": "repeats",
        "display_name": "Tandem / perfect repeats",
        "primary_domain": "advisory_scan",
        "order": 6,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Advisory finding.",
    },
    {
        "check_id": "gc_extremes",
        "display_name": "Local GC extremes (sliding window)",
        "primary_domain": "advisory_scan",
        "order": 7,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Local synthesis-hostile GC extremes, NOT the global GC optimization target.",
    },
    {
        "check_id": "splice_sites",
        "display_name": "Cryptic splice-site-like motifs",
        "primary_domain": "advisory_scan",
        "order": 8,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": (
            "Advisory pattern-match finding, low severity by default — not a validated "
            "splicing prediction."
        ),
    },
    {
        "check_id": "dinucleotides",
        "display_name": "CpG / TpA dinucleotide density",
        "primary_domain": "advisory_scan",
        "order": 9,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Advisory finding.",
    },
    {
        "check_id": "rare_codon_runs",
        "display_name": "Rare-codon runs",
        "primary_domain": "advisory_scan",
        "order": 10,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "advisory_only",
            "pipeline_default": "advisory_only",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": "Advisory finding, not a validated translation-rate measurement.",
    },
    # --- Domain D: Assembly Review -------------------------------------------
    {
        "check_id": "restriction_sites",
        "display_name": "Restriction Site Check (Type IIS)",
        "primary_domain": "assembly_review",
        "order": 11,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "not_checked",
            "pipeline_default": "hard_fail_raise_if_unfixable",
            "benchmark_scoring": "component_of_assembly_pass_soft",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": None}},
        "claim_boundary": (
            "Same biological check enforced at TWO different strengths: production pipeline "
            "raises and halts; benchmark scoring records a soft boolean."
        ),
    },
    {
        "check_id": "moclo_overhang",
        "display_name": "MoClo Level-0 overhang validity",
        "primary_domain": "assembly_review",
        "order": 12,
        "default_enabled": False,
        "enforcement_per_path": {
            "bare_optimizer": "not_checked",
            "pipeline_default": "opt_in_non_gating",
            "benchmark_scoring": "not_included",
        },
        "presentation": {"web": {"visible": True, "default_status_when_not_executed": "NOT_RUN"}},
        "claim_boundary": (
            "MUST be described as a non-gating, opt-in advisory check. It is not run by "
            "default, and even when run it never blocks a result — unlike restriction_sites, "
            "which can halt the pipeline."
        ),
    },
    # --- Domain A: Sequence Integrity (hard, gates pipeline via raise) ------
    {
        "check_id": "aa_identity",
        "display_name": "Amino-acid identity",
        "primary_domain": "integrity",
        "order": 13,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "not_checked",
            "pipeline_default": "hard_fail_raise",
            "benchmark_scoring": "component_of_biological_pass_soft",
        },
        "presentation": {"web": {"visible": False, "default_status_when_not_executed": None}},
        "claim_boundary": "Hard structural validity check, not a biological-risk prediction.",
    },
    {
        "check_id": "internal_stop",
        "display_name": "Internal stop codons",
        "primary_domain": "integrity",
        "order": 14,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "not_checked",
            "pipeline_default": "hard_fail_raise_via_final_validation",
            "benchmark_scoring": "component_of_biological_pass_soft",
        },
        "presentation": {"web": {"visible": False, "default_status_when_not_executed": None}},
        "claim_boundary": "Hard structural validity check.",
    },
    {
        "check_id": "invalid_codon",
        "display_name": "Invalid / partial codons",
        "primary_domain": "integrity",
        "order": 15,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "not_checked",
            "pipeline_default": "hard_fail_raise_via_final_validation",
            "benchmark_scoring": "component_of_biological_pass_soft",
        },
        "presentation": {"web": {"visible": False, "default_status_when_not_executed": None}},
        "claim_boundary": "Hard structural validity check.",
    },
    {
        "check_id": "frame_length",
        "display_name": "Reading-frame / length consistency",
        "primary_domain": "integrity",
        "order": 16,
        "default_enabled": True,
        "enforcement_per_path": {
            "bare_optimizer": "not_checked",
            "pipeline_default": "implicit_validator_stage",
            "benchmark_scoring": "component_of_biological_pass_soft",
        },
        "presentation": {"web": {"visible": False, "default_status_when_not_executed": None}},
        "claim_boundary": "Hard structural validity check.",
    },
    {
        "check_id": "custom_forbidden_motifs",
        "display_name": "User-configured forbidden motifs",
        "primary_domain": "configured_constraint",
        "order": 17,
        "default_enabled": False,
        "enforcement_per_path": {
            "bare_optimizer": "not_wired",
            "pipeline_default": "not_wired",
            "benchmark_scoring": "not_wired",
        },
        "presentation": {"web": {"visible": False, "default_status_when_not_executed": None}},
        "claim_boundary": (
            "MUST NOT be described as an active product feature. Implemented in a "
            "legacy/archived code path, not connected to the current optimize()/pipeline()/"
            "CLI/API surface at the pinned commit."
        ),
    },
)

REVIEW_CONTRACT_CHECK_IDS: tuple[str, ...] = (
    "polya",
    "are",
    "at_runs",
    "homopolymers",
    "repeats",
    "gc_extremes",
    "splice_sites",
    "dinucleotides",
    "rare_codon_runs",
    "restriction_sites",
    "moclo_overhang",
)

PUBLIC_VALIDATION_BADGE_IDS: tuple[str, ...] = ("global_gc_range",) + REVIEW_CONTRACT_CHECK_IDS

LEGACY_VALIDATION_FIELD_MAP: dict[str, str] = {
    "polya": "polya",
    "gc": "global_gc_range",
    # validation.moclo carries the Type IIS restriction-site result, NOT
    # moclo_overhang.
    "moclo": "restriction_sites",
}

_CHECKS_BY_ID: dict[str, dict[str, Any]] = {check["check_id"]: check for check in VALIDATION_CHECKS}


def get_check(check_id: str) -> dict[str, Any]:
    """Return the registry entry for `check_id`.

    Raises:
        KeyError: If `check_id` is not in the registry.
    """
    return _CHECKS_BY_ID[check_id]


def public_badge_checks() -> list[dict[str, Any]]:
    """Return the 12 public-badge checks, sorted by registry `order`."""
    return sorted(
        (check for check in VALIDATION_CHECKS if check["check_id"] in PUBLIC_VALIDATION_BADGE_IDS),
        key=lambda check: check["order"],
    )
