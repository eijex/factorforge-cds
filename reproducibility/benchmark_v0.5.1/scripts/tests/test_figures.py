import json
import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "figures"))
from make_benchmark_figures import load_summary, make_figure2, make_figure3, make_table3


MINIMAL_SUMMARY = {
    "factorforge_version": "3.2.0",
    "registry_version": "3.2.0",
    "dataset_n": 10,
    "random_seed": 320,
    "methods": {
        m: {
            "n_ok": 10, "n_error": 0, "error_rate": 0.0,
            "pass_rate_biological": 0.9, "pass_rate_assembly": 0.8,
            "pass_rate_multi_constraint": 0.6 if "assembly" in m else 0.3,
            "gc_in_range_rate": 0.95 if "gc" in m else 0.5,
            "mean_cai": 0.91 if "assembly" in m else 0.75,
        }
        for m in [
            "random_synonymous", "greedy_cai", "native_reference",
            "factorforge_balanced", "factorforge_gc_target",
            "factorforge_high_cai", "factorforge_assembly_friendly",
        ]
    },
}


def test_make_figure2_produces_files(tmp_path):
    make_figure2(MINIMAL_SUMMARY, tmp_path)
    assert (tmp_path / "figures/figure2_multiconstraint_pass_rate.png").exists()
    assert (tmp_path / "figures/figure2_multiconstraint_pass_rate.svg").exists()
    assert (tmp_path / "figures/figure2_multiconstraint_pass_rate.png").stat().st_size > 0


def test_make_figure3_produces_files(tmp_path):
    make_figure3(MINIMAL_SUMMARY, tmp_path)
    assert (tmp_path / "figures/figure3_benchmark_tradeoff_heatmap.png").exists()
    assert (tmp_path / "figures/figure3_benchmark_tradeoff_heatmap.svg").exists()


def test_make_table3_content(tmp_path):
    make_table3(MINIMAL_SUMMARY, tmp_path)
    table = (tmp_path / "tables/table3_benchmark_summary.md").read_text()
    assert "Multi-constraint pass" in table
    assert "FF assembly-friendly" in table
    assert "0.6000" in table  # assembly_friendly pass rate


def test_load_summary_real_file():
    real_path = Path("benchmarks/results/v3.2.0/benchmark_summary.json")
    if not real_path.exists():
        pytest.skip("Real benchmark_summary.json not available")
    summary = load_summary(real_path)
    assert summary["factorforge_version"] == "3.2.0"
    assert summary["dataset_n"] == 49257
    assert summary["random_seed"] == 320
    assert "factorforge_assembly_friendly" in summary["methods"]
