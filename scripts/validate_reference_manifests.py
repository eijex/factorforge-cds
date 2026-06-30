#!/usr/bin/env python3
"""Validate FactorForge codon-reference activation policy manifests.

This guard separates reference-table file presence from product activation.
It checks schema shape, file checksums, and default/candidate drift against the
current active codon-reference registry artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

try:  # Optional at runtime, required in the test environment.
    import jsonschema
except ImportError:  # pragma: no cover - exercised only in minimal envs.
    jsonschema = None  # type: ignore[assignment]

try:
    import yaml
except ImportError:  # pragma: no cover - exercised only in minimal envs.
    yaml = None  # type: ignore[assignment]

ROOT = Path(__file__).resolve().parents[1]
POLICY_MANIFEST_PATH = ROOT / "data" / "reference" / "reference_policy_manifest.json"
POLICY_SCHEMA_PATH = ROOT / "schemas" / "reference_policy_manifest.schema.json"
ACTIVE_REFERENCE_PATH = ROOT / "data" / "reference" / "active_codon_reference.json"
REGISTRY_PATH = ROOT / "src" / "factorforge" / "registry" / "current_parameter_registry.yaml"

ENABLED_STATUSES = {"enabled"}
NON_PRODUCTION_TIERS = {
    "historical_reproduction",
    "research_comparator",
    "experimental_candidate",
    "packaged_asset",
    "custom_user_provided",
}
TIER_ALLOWED_STATUSES = {
    "production_enabled": {"enabled"},
    "historical_reproduction": {"historical_only"},
    "research_comparator": {"research_only"},
    "experimental_candidate": {"disabled"},
    "packaged_asset": {"packaged_only", "disabled"},
    "custom_user_provided": {"disabled"},
}
REQUIRED_PRINCIPLES = {
    "reference file presence != product activation",
    "research comparator != public default",
    "disabled UI != backend permission",
    "manuscript framework != empirical validation claim",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_schema(manifest: dict[str, Any], schema: dict[str, Any]) -> None:
    if jsonschema is None:
        raise RuntimeError("jsonschema is required to validate reference policy manifests")
    jsonschema.validate(manifest, schema)


def load_registry() -> dict[str, Any]:
    if yaml is None:
        raise RuntimeError("PyYAML is required to validate registry/default drift")
    return yaml.safe_load(REGISTRY_PATH.read_text(encoding="utf-8"))


def reference_by_id(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    refs = manifest["references"]
    by_id = {entry["reference_id"]: entry for entry in refs}
    if len(by_id) != len(refs):
        raise AssertionError("reference_policy_manifest contains duplicate reference_id values")
    return by_id


def validate_policy_semantics(manifest: dict[str, Any], *, root: Path = ROOT) -> None:
    principles = set(manifest.get("activation_principles", []))
    missing = REQUIRED_PRINCIPLES - principles
    if missing:
        raise AssertionError(f"Missing activation principles: {sorted(missing)}")

    by_id = reference_by_id(manifest)
    active_id = manifest["active_default_reference_id"]
    if active_id not in by_id:
        raise AssertionError(f"active_default_reference_id {active_id!r} has no reference entry")

    enabled = [entry for entry in by_id.values() if entry["activation_status"] in ENABLED_STATUSES]
    if [entry["reference_id"] for entry in enabled] != [active_id]:
        raise AssertionError(
            "Exactly one enabled reference must exist and it must match active_default_reference_id; "
            f"enabled={[entry['reference_id'] for entry in enabled]!r}, active={active_id!r}"
        )

    for entry in by_id.values():
        tier = entry["tier"]
        status = entry["activation_status"]
        allowed = TIER_ALLOWED_STATUSES[tier]
        if status not in allowed:
            raise AssertionError(
                f"Reference {entry['reference_id']} has activation_status={status!r} "
                f"incompatible with tier={tier!r}; allowed={sorted(allowed)!r}"
            )
        if tier in NON_PRODUCTION_TIERS and status == "enabled":
            raise AssertionError(
                f"Reference {entry['reference_id']} is non-production tier {tier!r} but marked enabled"
            )

        codon_path = root / entry["codon_table_path"]
        if not codon_path.exists():
            raise AssertionError(f"Codon table path does not exist for {entry['reference_id']}: {codon_path}")
        actual = sha256_file(codon_path)
        expected = entry["checksum_sha256"]
        if actual != expected:
            raise AssertionError(
                f"Checksum mismatch for {entry['reference_id']}:\n  manifest: {expected}\n  actual:   {actual}"
            )


def validate_active_default_alignment(manifest: dict[str, Any]) -> None:
    active = load_json(ACTIVE_REFERENCE_PATH)
    registry = load_registry()
    policy_active_id = manifest["active_default_reference_id"]
    active_file_id = active["active_codon_table_id"]
    registry_active = registry["codon_reference"]["active"]

    if policy_active_id != active_file_id:
        raise AssertionError(
            f"Policy active_default_reference_id {policy_active_id!r} does not match "
            f"active_codon_reference.json {active_file_id!r}"
        )
    if policy_active_id != registry_active["id"]:
        raise AssertionError(
            f"Policy active_default_reference_id {policy_active_id!r} does not match "
            f"registry codon_reference.active.id {registry_active['id']!r}"
        )

    by_id = reference_by_id(manifest)
    policy_active = by_id[policy_active_id]
    if policy_active["codon_table_path"] != active["file"]:
        raise AssertionError(
            f"Policy active file {policy_active['codon_table_path']!r} does not match "
            f"active_codon_reference.json file {active['file']!r}"
        )
    if policy_active["checksum_sha256"] != registry_active["sha256"]:
        raise AssertionError(
            f"Policy active checksum {policy_active['checksum_sha256']!r} does not match "
            f"registry active sha256 {registry_active['sha256']!r}"
        )

    candidate = active.get("candidate", {})
    candidate_id = candidate.get("codon_table_id")
    if candidate_id:
        if candidate_id not in by_id:
            raise AssertionError(f"active_codon_reference candidate {candidate_id!r} missing from policy manifest")
        candidate_entry = by_id[candidate_id]
        if candidate_entry["activation_status"] == "enabled":
            raise AssertionError(f"Candidate {candidate_id!r} must not be production enabled")


def validate_all(
    manifest_path: Path = POLICY_MANIFEST_PATH,
    schema_path: Path = POLICY_SCHEMA_PATH,
) -> dict[str, Any]:
    manifest = load_json(manifest_path)
    schema = load_json(schema_path)
    validate_schema(manifest, schema)
    validate_policy_semantics(manifest)
    validate_active_default_alignment(manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=POLICY_MANIFEST_PATH)
    parser.add_argument("--schema", type=Path, default=POLICY_SCHEMA_PATH)
    args = parser.parse_args()

    manifest = validate_all(args.manifest, args.schema)
    print(
        "Reference policy manifest OK: "
        f"{len(manifest['references'])} references, active={manifest['active_default_reference_id']}"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
