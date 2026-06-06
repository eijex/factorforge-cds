"""A-stage guardrail: registry/spec integrity, not metric correctness."""
from __future__ import annotations
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "src" / "factorforge" / "registry" / "current_parameter_registry.yaml"
SPEC = ROOT / "benchmarks" / "benchmark_spec.yaml"


def _load(p): return yaml.safe_load(p.read_text(encoding="utf-8"))


def test_files_exist():
    assert REGISTRY.exists() and SPEC.exists()


def _iter_param_nodes(node, prefix=""):
    """Recursively yield (dotted_path, leaf_node) for all parameter leaf nodes."""
    if isinstance(node, dict) and "claim_level" in node:
        yield prefix, node
        return
    if isinstance(node, dict):
        for k, v in node.items():
            child_prefix = f"{prefix}.{k}" if prefix else k
            yield from _iter_param_nodes(v, child_prefix)


def test_every_param_has_required_fields():
    reg = _load(REGISTRY)
    required = {"claim_level", "evidence_status", "release_status", "rationale",
                "visibility", "permission", "provenance"}
    for path, body in _iter_param_nodes(reg["parameters"]):
        for field in required:
            assert field in body, f"{path} missing required field: {field}"


def test_claim_levels_valid():
    reg = _load(REGISTRY)
    allowed_claim = {"hard_requirement", "optimization_target", "warning_threshold",
                     "experimental_setting", "documentation_only"}
    for path, body in _iter_param_nodes(reg["parameters"]):
        assert body.get("claim_level") in allowed_claim, f"{path} bad claim_level"


def test_hard_requirements_have_hard_gate_evidence():
    reg = _load(REGISTRY)
    for path, body in _iter_param_nodes(reg["parameters"]):
        if body.get("claim_level") == "hard_requirement":
            assert body.get("evidence_status") == "hard_gate", (
                f"{path} hard_requirement must have evidence_status=hard_gate"
            )


def _resolve_ref(reg, dotted: str):
    node = reg
    for part in dotted.split("."):
        node = node[part]
    return node


def test_all_registry_refs_resolve():
    reg = _load(REGISTRY)
    spec = _load(SPEC)
    refs = []
    for layer in spec["hard_gates"].values():
        for item in layer.values():
            refs.append(item["registry_ref"])
    for item in spec["soft_targets"].values():
        refs.append(item["registry_ref"])
    for ref in refs:
        assert _resolve_ref(reg, ref) is not None, f"unresolved {ref}"


def test_hard_gates_split_into_two_classes():
    spec = _load(SPEC)
    assert set(spec["hard_gates"].keys()) == {"biological_validity", "assembly_validity"}
