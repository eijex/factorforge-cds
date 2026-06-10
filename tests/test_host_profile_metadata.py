"""Registry-backed host metadata guardrails for Open Bio Design Packages."""

import json
from pathlib import Path

import jsonschema

from factorforge.registry.registry_loader import load_registry, resolve_ref

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = json.loads(
    (ROOT / "src/factorforge/schemas/design_package.schema.json").read_text(encoding="utf-8")
)
REGISTRY = load_registry()
HOSTS = resolve_ref(REGISTRY, "parameters.host_profiles")


def test_nbenthamiana_metadata_is_stable_and_taxonomically_correct() -> None:
    host = HOSTS["nbenthamiana"]
    assert host["ncbi_taxonomy_id"] == 4100
    assert host["ncbi_taxonomy_id"] != 4197
    assert host["status"] == "stable"


def test_by2_remains_experimental() -> None:
    host = HOSTS["by2"]
    assert host["status"] == "experimental"
    assert host["status"] != "stable"


def test_wolffia_is_not_promoted_to_stable_host_profile() -> None:
    wolffia_profiles = {
        name: metadata for name, metadata in HOSTS.items() if "wolffia" in name.lower()
    }
    assert all(metadata["status"] in {"draft", "future"} for metadata in wolffia_profiles.values())


def test_schema_host_profiles_match_registry() -> None:
    alternatives = SCHEMA["$defs"]["hostProfile"]["oneOf"]
    schema_hosts = {
        item["properties"]["id"]["const"]: {
            key: value["const"] for key, value in item["properties"].items()
        }
        for item in alternatives
    }
    for host_id, schema_host in schema_hosts.items():
        registry_host = HOSTS[host_id]
        for field in (
            "display_name",
            "scientific_name",
            "ncbi_taxonomy_id",
            "status",
        ):
            assert schema_host[field] == registry_host[field]


def test_host_profile_roadmap_remains_deferred() -> None:
    roadmap = ROOT / "docs/host_profiles/HOST_PROFILE_REGISTRY_ROADMAP.md"
    assert not roadmap.exists(), "Host profile roadmap is deferred until Job 070 is Done"


def test_design_package_accepts_registry_host_metadata(design_package: dict) -> None:
    jsonschema.validate(design_package, SCHEMA)
