"""Load and resolve the factorforge parameter registry (package source of truth)."""
from __future__ import annotations
from pathlib import Path
import yaml

REGISTRY_PATH = Path(__file__).resolve().parent / "current_parameter_registry.yaml"


def load_registry() -> dict:
    return yaml.safe_load(REGISTRY_PATH.read_text(encoding="utf-8"))


def resolve_ref(registry: dict, dotted: str):
    """Resolve a dotted path (e.g. 'parameters.optimization.cai_target') into a value."""
    node = registry
    for part in dotted.split("."):
        node = node[part]
    return node
