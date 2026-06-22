"""Codon table manifest provenance guardrails."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "data" / "reference" / "codon_table_manifest.json"
SCHEMA_PATH = ROOT / "schemas" / "codon_table_manifest.schema.json"
CODON_TABLE_PATH = ROOT / "src" / "factorforge" / "data" / "nbenthamiana_codons.json"
ARCHIVE_SCRIPT_PATH = ROOT / "archive" / "v3-ml-prototype" / "scripts" / "1_data_preparation" / "build_codon_table.py"


@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def schema():
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def test_manifest_and_schema_exist():
    assert MANIFEST_PATH.exists(), f"Missing: {MANIFEST_PATH}"
    assert SCHEMA_PATH.exists(), f"Missing: {SCHEMA_PATH}"


def test_manifest_required_fields(manifest):
    required = [
        "codon_table_id",
        "sha256",
        "source_status",
        "build_path_status",
        "authoritative_build_script_available",
    ]
    for field in required:
        assert field in manifest, f"Missing required field: {field}"


def test_manifest_authoritative_build_script_false(manifest):
    assert manifest["authoritative_build_script_available"] is False, (
        "Current legacy codon table must declare authoritative_build_script_available: false"
    )


def test_manifest_build_path_incomplete(manifest):
    assert manifest["build_path_status"] == "incomplete", (
        f"Expected build_path_status='incomplete', got {manifest['build_path_status']!r}"
    )


def test_manifest_source_status_legacy(manifest):
    assert manifest["source_status"] == "legacy_metadata_only", (
        f"Expected source_status='legacy_metadata_only', got {manifest['source_status']!r}"
    )


def test_manifest_sha256_matches_actual_file(manifest):
    assert CODON_TABLE_PATH.exists(), f"Codon table file not found: {CODON_TABLE_PATH}"
    actual = hashlib.sha256(CODON_TABLE_PATH.read_bytes()).hexdigest()
    assert manifest["sha256"] == actual, (
        f"SHA-256 mismatch:\n  manifest: {manifest['sha256']}\n  actual:   {actual}\n"
        "Update manifest sha256 if the codon table was intentionally changed."
    )


def test_archive_script_only_in_not_authoritative(manifest):
    """Archive ML prototype script must NOT be cited as the authoritative builder."""
    archive_rel = str(ARCHIVE_SCRIPT_PATH.relative_to(ROOT)).replace("\\", "/")

    # Must not appear as a top-level authoritative field
    for key in ("authoritative_build_script", "build_script", "source_script"):
        assert key not in manifest, (
            f"Field '{key}' must not exist in manifest — "
            "archive script must not be cited as authoritative."
        )

    # If mentioned at all, must only be in not_authoritative_build_scripts
    not_auth = manifest.get("not_authoritative_build_scripts", [])
    not_auth_paths = {e.get("path", "").replace("\\", "/") for e in not_auth}
    assert archive_rel in not_auth_paths, (
        f"Archive script {archive_rel!r} must be listed in not_authoritative_build_scripts"
    )


def test_schema_required_fields_present(schema):
    required = schema.get("required", [])
    for field in ["codon_table_id", "sha256", "source_status", "build_path_status",
                  "authoritative_build_script_available"]:
        assert field in required, f"Schema missing required field: {field}"


def test_manifest_validates_against_schema(manifest, schema):
    pytest.importorskip("jsonschema", reason="jsonschema not installed — skipping live validation")
    import jsonschema
    jsonschema.validate(manifest, schema)


def test_ncbi_taxid_nbenthamiana(manifest):
    assert manifest["ncbi_taxid"] == 4100, (
        f"N. benthamiana taxid must be 4100, got {manifest['ncbi_taxid']}"
    )
