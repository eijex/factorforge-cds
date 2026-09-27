"""Single source of truth for FactorForge product and engine versions."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files
from typing import Any


@lru_cache(maxsize=1)
def version_manifest() -> dict[str, Any]:
    manifest_path = files("factorforge.registry").joinpath("version_manifest.json")
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def product_version() -> str:
    return str(version_manifest()["product"]["version"])


def engine_metadata(engine: str) -> dict[str, Any]:
    return dict(version_manifest()["engines"][engine])


def engine_version(engine: str) -> str:
    return str(engine_metadata(engine)["version"])


def engine_generation(engine: str) -> int:
    return int(engine_metadata(engine).get("generation", 1))


def engine_status(engine: str) -> str:
    return str(engine_metadata(engine).get("status", "unknown"))


def engine_runtime_state(engine: str) -> str:
    return str(engine_metadata(engine).get("runtime_state", "stable"))


def public_version_metadata() -> dict[str, Any]:
    manifest = version_manifest()
    return {
        "product": dict(manifest["product"]),
        "engines": {name: dict(metadata) for name, metadata in manifest["engines"].items()},
    }
