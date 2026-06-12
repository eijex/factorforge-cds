import sys
import json
from pathlib import Path
import yaml
import pytest
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[4]
SPEC_PATH = REPO_ROOT / "benchmarks" / "ablation" / "ablation_spec.yaml"

sys.path.insert(0, str(REPO_ROOT))

from benchmarks.ablation.conditions.cai_gc import ablation_cai_gc_cds
from benchmarks.ablation.conditions.cai_type_iis import ablation_cai_type_iis_cds
from benchmarks.ablation.conditions.cai_gc_type_iis import ablation_cai_gc_type_iis_cds
from benchmarks.scoring import score_cds
from benchmarks.config import load_benchmark_config

_CFG = load_benchmark_config()
_TEST_PROTEIN = "MARNKV"  # 6 AAs
FORBIDDEN_PATTERNS = ["GGTCTC", "GAGACC", "GAAGAC", "GTCTTC", "CGTCTC", "GAGACG"]


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


def test_cai_gc_returns_correct_length():
    cds = ablation_cai_gc_cds(_TEST_PROTEIN)
    assert len(cds) == len(_TEST_PROTEIN) * 3


def test_cai_gc_aa_identity():
    cds = ablation_cai_gc_cds(_TEST_PROTEIN)
    scores = score_cds("ablation_cai_gc", "ablation", "t0", _TEST_PROTEIN, cds, _CFG, 0.0)
    assert scores["aa_identity"] == 1.0


def test_cai_gc_is_deterministic():
    cds1 = ablation_cai_gc_cds(_TEST_PROTEIN)
    cds2 = ablation_cai_gc_cds(_TEST_PROTEIN)
    assert cds1 == cds2


def test_cai_type_iis_returns_correct_length():
    cds = ablation_cai_type_iis_cds(_TEST_PROTEIN)
    assert len(cds) == len(_TEST_PROTEIN) * 3


def test_cai_type_iis_aa_identity():
    cds = ablation_cai_type_iis_cds(_TEST_PROTEIN)
    scores = score_cds("ablation_cai_type_iis", "ablation", "t0", _TEST_PROTEIN, cds, _CFG, 0.0)
    assert scores["aa_identity"] == 1.0


def test_cai_type_iis_no_forbidden_sites_short_protein():
    cds = ablation_cai_type_iis_cds(_TEST_PROTEIN, seed=320, max_attempts=50)
    for pat in FORBIDDEN_PATTERNS:
        assert pat not in cds, f"Found forbidden pattern {pat} in {cds}"


def test_cai_type_iis_seed_reproducibility():
    cds1 = ablation_cai_type_iis_cds(_TEST_PROTEIN, seed=42)
    cds2 = ablation_cai_type_iis_cds(_TEST_PROTEIN, seed=42)
    assert cds1 == cds2


def test_cai_gc_type_iis_returns_correct_length():
    cds = ablation_cai_gc_type_iis_cds(_TEST_PROTEIN)
    assert len(cds) == len(_TEST_PROTEIN) * 3


def test_cai_gc_type_iis_aa_identity():
    cds = ablation_cai_gc_type_iis_cds(_TEST_PROTEIN)
    scores = score_cds("ablation_cai_gc_type_iis", "ablation", "t0", _TEST_PROTEIN, cds, _CFG, 0.0)
    assert scores["aa_identity"] == 1.0


def test_cai_gc_type_iis_no_forbidden_sites_short_protein():
    cds = ablation_cai_gc_type_iis_cds(_TEST_PROTEIN, seed=320, max_attempts=50)
    for pat in FORBIDDEN_PATTERNS:
        assert pat not in cds, f"Found forbidden pattern {pat} in {cds}"


def test_cai_gc_type_iis_seed_reproducibility():
    cds1 = ablation_cai_gc_type_iis_cds(_TEST_PROTEIN, seed=77)
    cds2 = ablation_cai_gc_type_iis_cds(_TEST_PROTEIN, seed=77)
    assert cds1 == cds2


def test_cai_gc_type_iis_uses_codon_weights():
    # Verify the module-level _WEIGHTS attribute exists (CAI component required).
    # If L4 used GC-only weights, it would not need _WEIGHTS — this guards against regression.
    import benchmarks.ablation.conditions.cai_gc_type_iis as m
    assert hasattr(m, "_WEIGHTS"), "L4 module must import and use _WEIGHTS for CAI weighting"


# --- Fixture mini CSV ---
def _make_mini_formal_csv(tmp_dir) -> Path:
    """Create a minimal benchmark_results.csv fixture with 3 methods × 3 sequences."""
    rows = []
    for sid in ["seq001", "seq002", "seq003"]:
        for rep in [1, 2, 3]:
            rows.append({
                "method": "random_synonymous", "method_type": "baseline",
                "sequence_id": sid, "replicate": rep, "seed": 320,
                "aa_identity": 1.0, "internal_stop_count": 0, "invalid_codon_count": 0,
                "length_multiple_of_three": True, "cai": 0.50,
                "gc_percent": 58.0, "gc_in_target_range": True,
                "forbidden_type_iis_site_count": 0,
                "biological_pass": True, "assembly_pass": True,
                "multi_constraint_pass": True, "runtime_seconds": 0.01,
                "status": "ok", "error_type": None, "error_message": None,
            })
        rows.append({
            "method": "greedy_cai", "method_type": "baseline",
            "sequence_id": sid, "replicate": 1, "seed": "",
            "aa_identity": 1.0, "internal_stop_count": 0, "invalid_codon_count": 0,
            "length_multiple_of_three": True, "cai": 0.75,
            "gc_percent": 62.0, "gc_in_target_range": True,
            "forbidden_type_iis_site_count": 1,
            "biological_pass": True, "assembly_pass": False,
            "multi_constraint_pass": False, "runtime_seconds": 0.01,
            "status": "ok", "error_type": None, "error_message": None,
        })
        rows.append({
            "method": "factorforge_assembly_friendly", "method_type": "optimizer",
            "sequence_id": sid, "replicate": 1, "seed": 320,
            "aa_identity": 1.0, "internal_stop_count": 0, "invalid_codon_count": 0,
            "length_multiple_of_three": True, "cai": 0.70,
            "gc_percent": 59.0, "gc_in_target_range": True,
            "forbidden_type_iis_site_count": 0,
            "biological_pass": True, "assembly_pass": True,
            "multi_constraint_pass": True, "runtime_seconds": 0.05,
            "status": "ok", "error_type": None, "error_message": None,
        })
    p = Path(tmp_dir) / "benchmark_results.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_extract_l0_from_formal_csv(tmp_path):
    from benchmarks.ablation.run_ablation import extract_existing_layers
    csv_path = _make_mini_formal_csv(tmp_path)
    df = pd.read_csv(csv_path)
    l0 = extract_existing_layers(df, "L0", "random_synonymous",
                                  enabled_constraints={"cai": False, "gc_target": False, "type_iis_clean": False})
    assert len(l0) == 9  # 3 seqs × 3 replicates
    assert (l0["ablation_layer"] == "L0").all()
    assert (l0["ablation_source"] == "existing_csv").all()
    assert "enabled_constraints_json" in l0.columns


def test_extract_l1_from_formal_csv(tmp_path):
    from benchmarks.ablation.run_ablation import extract_existing_layers
    csv_path = _make_mini_formal_csv(tmp_path)
    df = pd.read_csv(csv_path)
    l1 = extract_existing_layers(df, "L1", "greedy_cai",
                                  enabled_constraints={"cai": True, "gc_target": False, "type_iis_clean": False})
    assert len(l1) == 3  # 3 seqs
    assert (l1["ablation_layer"] == "L1").all()
    assert "enabled_constraints_json" in l1.columns


def test_limit_filters_existing_layers_to_same_sequence_ids(tmp_path):
    from benchmarks.ablation.run_ablation import extract_existing_layers
    csv_path = _make_mini_formal_csv(tmp_path)
    full_df = pd.read_csv(csv_path)
    subset_ids = {"seq001", "seq002"}
    filtered_df = full_df[full_df["sequence_id"].isin(subset_ids)].copy()
    l0 = extract_existing_layers(filtered_df, "L0", "random_synonymous",
                                  enabled_constraints={"cai": False, "gc_target": False, "type_iis_clean": False})
    assert set(l0["sequence_id"].unique()) == subset_ids
    assert len(l0) == 6  # 2 seqs × 3 replicates


def test_ablation_summary_has_required_fields(tmp_path):
    from benchmarks.ablation.run_ablation import build_ablation_summary
    spec = yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))
    rows = []
    for layer in ["L0", "L1", "L5"]:
        for sid in ["seq001", "seq002"]:
            rows.append({
                "ablation_layer": layer, "ablation_source": "existing_csv",
                "sequence_id": sid, "multi_constraint_pass": True,
                "biological_pass": True, "assembly_pass": True,
                "gc_in_target_range": True, "cai": 0.70, "gc_percent": 60.0,
                "forbidden_type_iis_site_count": 0, "status": "ok",
            })
    df = pd.DataFrame(rows)
    summary = build_ablation_summary(df, spec, "abc123", "def456", "ghi789", "jkl012")
    assert summary["analysis_type"] == "constraint_ablation"
    assert "source_formal_run_id" in summary
    assert "layers" in summary
    for layer in ["L0", "L1", "L5"]:
        assert layer in summary["layers"]
