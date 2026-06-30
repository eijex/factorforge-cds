"""Reference policy manifest guardrails.

These tests enforce the Job 176/177 split between packaged codon-reference
assets and production activation.
"""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
VALIDATOR_PATH = ROOT / "scripts" / "validate_reference_manifests.py"
POLICY_MANIFEST_PATH = ROOT / "data" / "reference" / "reference_policy_manifest.json"
POLICY_SCHEMA_PATH = ROOT / "schemas" / "reference_policy_manifest.schema.json"

spec = importlib.util.spec_from_file_location("validate_reference_manifests", VALIDATOR_PATH)
validator = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(validator)


@pytest.fixture(scope="module")
def policy_manifest() -> dict:
    return json.loads(POLICY_MANIFEST_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def policy_schema() -> dict:
    return json.loads(POLICY_SCHEMA_PATH.read_text(encoding="utf-8"))


def test_reference_policy_manifest_and_schema_exist() -> None:
    assert POLICY_MANIFEST_PATH.exists(), f"Missing: {POLICY_MANIFEST_PATH}"
    assert POLICY_SCHEMA_PATH.exists(), f"Missing: {POLICY_SCHEMA_PATH}"
    assert VALIDATOR_PATH.exists(), f"Missing: {VALIDATOR_PATH}"


def test_reference_policy_manifest_validates_against_schema(policy_manifest, policy_schema) -> None:
    jsonschema = pytest.importorskip("jsonschema", reason="jsonschema not installed")
    jsonschema.validate(policy_manifest, policy_schema)


def test_reference_policy_manifest_passes_live_validator() -> None:
    manifest = validator.validate_all(POLICY_MANIFEST_PATH, POLICY_SCHEMA_PATH)
    assert manifest["active_default_reference_id"] == "nbenthamiana_legacy_kazusa_sgn_v101"


def test_required_activation_principles_are_present(policy_manifest) -> None:
    principles = set(policy_manifest["activation_principles"])
    for principle in validator.REQUIRED_PRINCIPLES:
        assert principle in principles


def test_reference_ids_are_unique(policy_manifest) -> None:
    refs = [entry["reference_id"] for entry in policy_manifest["references"]]
    assert len(refs) == len(set(refs))


def test_only_active_default_is_enabled(policy_manifest) -> None:
    active = policy_manifest["active_default_reference_id"]
    enabled = [
        entry["reference_id"]
        for entry in policy_manifest["references"]
        if entry["activation_status"] == "enabled"
    ]
    assert enabled == [active]


def test_non_production_references_are_not_enabled(policy_manifest) -> None:
    for entry in policy_manifest["references"]:
        if entry["tier"] != "production_enabled":
            assert entry["activation_status"] != "enabled", entry["reference_id"]


def test_candidate_reference_is_disabled(policy_manifest) -> None:
    refs = {entry["reference_id"]: entry for entry in policy_manifest["references"]}
    candidate = refs["nbenthamiana_nbev11_hc_v2"]
    assert candidate["tier"] == "experimental_candidate"
    assert candidate["activation_status"] == "disabled"
    assert "not wet-lab validation" in candidate["claim_boundary"].lower()


def test_research_comparators_are_research_only(policy_manifest) -> None:
    refs = {entry["reference_id"]: entry for entry in policy_manifest["references"]}
    for ref_id in ["nbenthamiana_nbev11_all_cds_v2", "nbenthamiana_qld183_v103"]:
        assert refs[ref_id]["tier"] == "research_comparator"
        assert refs[ref_id]["activation_status"] == "research_only"
        assert "not accepted" in refs[ref_id]["public_exposure"]["rest_api"]


def test_checksum_drift_is_detected(policy_manifest) -> None:
    mutated = copy.deepcopy(policy_manifest)
    mutated["references"][0]["checksum_sha256"] = "0" * 64
    with pytest.raises(AssertionError, match="Checksum mismatch"):
        validator.validate_policy_semantics(mutated)


def test_missing_active_reference_is_detected(policy_manifest) -> None:
    mutated = copy.deepcopy(policy_manifest)
    mutated["active_default_reference_id"] = "missing_reference"
    with pytest.raises(AssertionError, match="has no reference entry"):
        validator.validate_policy_semantics(mutated)


def test_research_comparator_cannot_be_enabled(policy_manifest) -> None:
    mutated = copy.deepcopy(policy_manifest)
    for entry in mutated["references"]:
        if entry["tier"] == "research_comparator":
            entry["activation_status"] = "enabled"
            break
    with pytest.raises(AssertionError, match="enabled reference|incompatible with tier"):
        validator.validate_policy_semantics(mutated)


def test_active_default_alignment_with_registry_and_active_reference(policy_manifest) -> None:
    validator.validate_active_default_alignment(policy_manifest)
