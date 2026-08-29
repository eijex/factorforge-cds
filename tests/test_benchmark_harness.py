import pytest
import os
import json
from factorforge.benchmark.models import BenchmarkRunConfig, BenchmarkTarget
from factorforge.benchmark.runner import BenchmarkRunner

def test_benchmark_runner_engineering_suite():
    # Setup configuration
    config = BenchmarkRunConfig(
        suite_name="Engineering Benchmark Suite Test",
        host="nbenthamiana",
        target_gc_min_percent=40.0,
        target_gc_max_percent=47.0,
        terminal_stop_policy="append",
        forbidden_type_iis=["BsaI"]
    )
    
    # Tiny target for testing
    targets = [
        BenchmarkTarget(target_id="TinyTest", protein_sequence="MVLSA")
    ]
    
    # We will test all 3 arms (Profile, DP stub, LM scaffold)
    runner = BenchmarkRunner(config)
    report = runner.run_suite(targets, ["profile", "dp", "lm"])
    
    # Assertions
    assert report.suite_name == "Engineering Benchmark Suite Test"
    assert len(report.results) == 3 # 1 target * 3 engines
    
    # Verify provenance data
    res_lm = next(r for r in report.results if r.engine_name == "lm")
    assert res_lm.provenance.input_type == "protein"
    assert res_lm.provenance.constraint_profile_id == "default_gc"
    assert res_lm.provenance.terminal_stop_policy == "append"
    
    # Verify raw counts pass rate
    assert "profile" in report.pass_counts
    assert "dp" in report.pass_counts
    assert "lm" in report.pass_counts
    assert report.pass_counts["profile"].endswith("/1")
    
    # Check that average metrics exist
    assert "lm" in report.average_cai
    assert "dp" in report.average_gc
