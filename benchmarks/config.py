"""Load benchmark spec; resolve registry_ref values via factorforge.registry.registry_loader.
benchmarks/config.py is a consumer only — it does NOT own or reload the registry."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import yaml

from factorforge.registry.registry_loader import load_registry, resolve_ref, REGISTRY_PATH

ROOT = Path(__file__).resolve().parents[1]


def _load_spec(p: Path) -> dict:
    return yaml.safe_load(p.read_text(encoding="utf-8"))


@dataclass(frozen=True)
class BenchmarkConfig:
    gc_min: float
    gc_max: float
    cai_target: float
    forbidden_type_iis: list[str]
    random_seed: int
    random_replicates: int
    registry_version: str
    registry_sha256: str
    spec_sha256: str


def load_benchmark_config() -> BenchmarkConfig:
    import hashlib
    spec_path = ROOT / "benchmarks" / "benchmark_spec.yaml"
    spec = _load_spec(spec_path)
    registry = load_registry()
    gc_node = resolve_ref(registry, spec["soft_targets"]["gc_range"]["registry_ref"])
    gc_min, gc_max = gc_node["value"]
    cai = resolve_ref(registry, spec["soft_targets"]["cai"]["registry_ref"])
    type_iis = resolve_ref(registry, spec["hard_gates"]["assembly_validity"]["forbidden_type_iis_sites"]["registry_ref"])
    return BenchmarkConfig(
        gc_min=float(gc_min), gc_max=float(gc_max),
        cai_target=float(cai["value"]),
        forbidden_type_iis=list(type_iis["value"]),
        random_seed=int(spec["random"]["seed"]),
        random_replicates=int(spec["random"]["replicates"]),
        registry_version=registry["registry"]["version"],
        registry_sha256=hashlib.sha256(REGISTRY_PATH.read_bytes()).hexdigest(),
        spec_sha256=hashlib.sha256(spec_path.read_bytes()).hexdigest(),
    )
