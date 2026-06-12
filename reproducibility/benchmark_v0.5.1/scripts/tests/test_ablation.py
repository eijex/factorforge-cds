from pathlib import Path
import yaml
import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]
SPEC_PATH = REPO_ROOT / "benchmarks" / "ablation" / "ablation_spec.yaml"


def test_ablation_spec_parses():
    spec = yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))
    assert "ablation" in spec
    assert "layers" in spec
    layers = spec["layers"]
    assert set(layers.keys()) == {"L0", "L1", "L2", "L3", "L4", "L5"}


def test_ablation_spec_l0_l1_l5_are_existing_csv():
    spec = yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))
    for layer in ["L0", "L1", "L5"]:
        assert spec["layers"][layer]["source"] == "existing_csv"


def test_ablation_spec_l2_l3_l4_are_new_run():
    spec = yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))
    for layer in ["L2", "L3", "L4"]:
        assert spec["layers"][layer]["source"] == "new_run"
        assert "condition_fn" in spec["layers"][layer]


def test_ablation_spec_all_layers_have_enabled_constraints():
    spec = yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))
    for layer_name, layer in spec["layers"].items():
        assert "enabled_constraints" in layer, f"{layer_name} missing enabled_constraints"
