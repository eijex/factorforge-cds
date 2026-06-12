import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "analysis"))
from failure_attribution import categorize_row, aggregate_method, compute_overlap_matrix
from common import FAILURE_CATEGORIES


def _row(gc_in_range=True, iis_count=0, bio_pass=True, mc_pass=True, status="ok"):
    return {
        "gc_in_target_range": gc_in_range,
        "forbidden_type_iis_site_count": iis_count,
        "biological_pass": bio_pass,
        "multi_constraint_pass": mc_pass,
        "status": status,
    }


def test_categorize_passed_all():
    assert categorize_row(_row()) == "passed_all"


def test_categorize_failed_gc_only():
    assert categorize_row(_row(gc_in_range=False, mc_pass=False)) == "failed_gc_only"


def test_categorize_failed_iis_only():
    assert categorize_row(_row(iis_count=1, mc_pass=False)) == "failed_iis_only"


def test_categorize_failed_gc_and_iis():
    assert categorize_row(_row(gc_in_range=False, iis_count=2, mc_pass=False)) == "failed_gc_and_iis"


def test_categorize_failed_bio():
    assert categorize_row(_row(bio_pass=False, mc_pass=False)) == "failed_bio"


def test_categorize_failed_other_bad_status():
    assert categorize_row(_row(status="error")) == "failed_other"


def test_aggregate_deterministic_method():
    df = pd.DataFrame([
        {"sequence_id": f"seq{i}", "replicate": 1,
         "gc_in_target_range": True, "forbidden_type_iis_site_count": 0,
         "biological_pass": True, "multi_constraint_pass": True, "status": "ok"}
        for i in range(8)
    ] + [
        {"sequence_id": f"fail{i}", "replicate": 1,
         "gc_in_target_range": False, "forbidden_type_iis_site_count": 0,
         "biological_pass": True, "multi_constraint_pass": False, "status": "ok"}
        for i in range(2)
    ])
    result = aggregate_method(df, "greedy_cai")
    assert result["passed_all"] == 8.0
    assert result["failed_gc_only"] == 2.0
    assert result["total_sequence_equivalents"] == 10.0
    assert abs(result["passed_all_rate"] - 0.8) < 1e-9
    total = sum(result[cat] for cat in FAILURE_CATEGORIES)
    assert abs(total - 10.0) < 1e-6


def test_aggregate_replicated_method_weighted():
    # 2 sequences × 3 replicates
    # seq0: rep1=pass, rep2=fail_gc, rep3=pass → contributes 2/3 pass, 1/3 fail_gc
    # seq1: rep1=fail_iis, rep2=fail_iis, rep3=fail_iis → contributes 1.0 fail_iis
    rows = [
        {"sequence_id": "seq0", "replicate": 1, "gc_in_target_range": True,
         "forbidden_type_iis_site_count": 0, "biological_pass": True,
         "multi_constraint_pass": True, "status": "ok"},
        {"sequence_id": "seq0", "replicate": 2, "gc_in_target_range": False,
         "forbidden_type_iis_site_count": 0, "biological_pass": True,
         "multi_constraint_pass": False, "status": "ok"},
        {"sequence_id": "seq0", "replicate": 3, "gc_in_target_range": True,
         "forbidden_type_iis_site_count": 0, "biological_pass": True,
         "multi_constraint_pass": True, "status": "ok"},
        {"sequence_id": "seq1", "replicate": 1, "gc_in_target_range": True,
         "forbidden_type_iis_site_count": 1, "biological_pass": True,
         "multi_constraint_pass": False, "status": "ok"},
        {"sequence_id": "seq1", "replicate": 2, "gc_in_target_range": True,
         "forbidden_type_iis_site_count": 1, "biological_pass": True,
         "multi_constraint_pass": False, "status": "ok"},
        {"sequence_id": "seq1", "replicate": 3, "gc_in_target_range": True,
         "forbidden_type_iis_site_count": 1, "biological_pass": True,
         "multi_constraint_pass": False, "status": "ok"},
    ]
    df = pd.DataFrame(rows)
    result = aggregate_method(df, "random_synonymous")
    assert abs(result["passed_all"] - 2/3) < 1e-6
    assert abs(result["failed_gc_only"] - 1/3) < 1e-6
    assert abs(result["failed_iis_only"] - 1.0) < 1e-6
    assert abs(result["total_sequence_equivalents"] - 2.0) < 1e-6
    total = sum(result[cat] for cat in FAILURE_CATEGORIES)
    assert abs(total - 2.0) < 1e-6


def test_overlap_matrix_sums_deterministic():
    df = pd.DataFrame([
        {"method": "greedy_cai", "sequence_id": f"s{i}", "replicate": 1,
         "gc_in_target_range": i % 2 == 0, "forbidden_type_iis_site_count": i % 3}
        for i in range(9)
    ])
    matrix = compute_overlap_matrix(df)
    row = matrix[matrix["method"] == "greedy_cai"].iloc[0]
    total = row["gc_fail_only"] + row["iis_fail_only"] + row["gc_and_iis_fail"] + row["neither_fail"]
    assert total == 9
    assert row["n_sequence_equivalents"] == 9.0


def test_overlap_matrix_sums_replicated():
    # 3 sequences x 3 replicates for random_synonymous
    rows = [
        {"method": "random_synonymous", "sequence_id": f"s{i}", "replicate": rep,
         "gc_in_target_range": True, "forbidden_type_iis_site_count": 0}
        for i in range(3) for rep in range(1, 4)
    ]
    df = pd.DataFrame(rows)
    matrix = compute_overlap_matrix(df)
    row = matrix[matrix["method"] == "random_synonymous"].iloc[0]
    total = row["gc_fail_only"] + row["iis_fail_only"] + row["gc_and_iis_fail"] + row["neither_fail"]
    assert abs(total - 3.0) < 1e-6  # 3 sequence-equivalents
    assert row["n_sequence_equivalents"] == 3.0
