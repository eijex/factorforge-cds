"""
Failure attribution from benchmark_results.csv.
Handles random_synonymous 3-replicate weighted aggregation.
Run from repo root:
  python reproducibility/benchmark_v0.5.1/scripts/analysis/failure_attribution.py
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from common import METHOD_ORDER, FAILURE_CATEGORIES, REPLICATED_METHOD


def categorize_row(row) -> str:
    """Assign failure category to one benchmark result row."""
    status = str(row.get("status", "ok")).strip()
    if status not in ("ok", ""):
        return "failed_other"
    bio_fail = not bool(row["biological_pass"])
    gc_fail = not bool(row["gc_in_target_range"])
    iis_fail = int(row["forbidden_type_iis_site_count"]) > 0
    if bio_fail:
        return "failed_bio"
    if bool(row["multi_constraint_pass"]):
        return "passed_all"
    if gc_fail and iis_fail:
        return "failed_gc_and_iis"
    if gc_fail:
        return "failed_gc_only"
    if iis_fail:
        return "failed_iis_only"
    return "failed_other"


def aggregate_method(df_method: pd.DataFrame, method: str) -> dict:
    """
    Aggregate failure categories for one method.
    For random_synonymous: replicate-weighted (weight=1/n_replicates per sequence).
    For deterministic methods: direct count.
    Returns counts summing to 49257 sequence-equivalents and corresponding rates.
    """
    df_method = df_method.copy()
    df_method["category"] = df_method.apply(categorize_row, axis=1)

    if method == REPLICATED_METHOD and "replicate" in df_method.columns:
        n_replicates = df_method["replicate"].nunique()
        if n_replicates > 1:
            # Each (sequence_id, replicate) pair contributes 1/n_replicates
            weighted = (
                df_method.groupby(["sequence_id", "category"])
                .size()
                .unstack(fill_value=0)
                .div(n_replicates)
            )
            counts = {cat: float(weighted[cat].sum()) if cat in weighted.columns else 0.0
                      for cat in FAILURE_CATEGORIES}
            total = float(df_method["sequence_id"].nunique())
            return {**counts,
                    **{f"{cat}_rate": counts[cat] / total for cat in FAILURE_CATEGORIES},
                    "total_sequence_equivalents": total,
                    "n_replicates": n_replicates,
                    "aggregation": "replicate_weighted"}

    counts = df_method["category"].value_counts().to_dict()
    counts = {cat: float(counts.get(cat, 0)) for cat in FAILURE_CATEGORIES}
    total = float(len(df_method))
    return {**counts,
            **{f"{cat}_rate": counts[cat] / total for cat in FAILURE_CATEGORIES},
            "total_sequence_equivalents": total,
            "n_replicates": 1,
            "aggregation": "direct"}


def compute_overlap_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    GC x Type IIS co-failure counts per method, replicate-weighted.
    random_synonymous: each replicate contributes 1/n_replicates weight,
    so totals sum to 49,257 sequence-equivalents (consistent with failure attribution).
    """
    rows = []
    for method, grp in df.groupby("method"):
        n_replicates = grp["replicate"].nunique() if "replicate" in grp.columns else 1
        is_replicated = (method == REPLICATED_METHOD and n_replicates > 1)
        gc_fail = ~grp["gc_in_target_range"].astype(bool)
        iis_fail = grp["forbidden_type_iis_site_count"].astype(int) > 0
        if is_replicated:
            w = 1.0 / n_replicates
            rows.append({
                "method": method,
                "n_sequence_equivalents": float(grp["sequence_id"].nunique()),
                "gc_fail_only": float((gc_fail & ~iis_fail).sum()) * w,
                "iis_fail_only": float((~gc_fail & iis_fail).sum()) * w,
                "gc_and_iis_fail": float((gc_fail & iis_fail).sum()) * w,
                "neither_fail": float((~gc_fail & ~iis_fail).sum()) * w,
            })
        else:
            rows.append({
                "method": method,
                "n_sequence_equivalents": float(len(grp)),
                "gc_fail_only": int((gc_fail & ~iis_fail).sum()),
                "iis_fail_only": int((~gc_fail & iis_fail).sum()),
                "gc_and_iis_fail": int((gc_fail & iis_fail).sum()),
                "neither_fail": int((~gc_fail & ~iis_fail).sum()),
            })
    return pd.DataFrame(rows)


def run(csv_path: Path, out_dir: Path) -> None:
    print(f"Loading {csv_path} ...")
    df = pd.read_csv(csv_path)

    # Input audit
    has_status = "status" in df.columns
    audit = {
        "total_rows": len(df),
        "methods": sorted(df["method"].unique().tolist()),
        "has_status_column": has_status,
        "row_counts_per_method": {k: int(v) for k, v in df.groupby("method").size().items()},
        "unique_sequences": int(df["sequence_id"].nunique()),
        "dataset_n": int(df["sequence_id"].nunique()),
    }
    (out_dir / "data").mkdir(parents=True, exist_ok=True)
    (out_dir / "data/input_audit.json").write_text(json.dumps(audit, indent=2))
    print(f"  Input audit: {len(audit['methods'])} methods, {audit['unique_sequences']} sequences")

    # Failure attribution
    results = {}
    for method in METHOD_ORDER:
        if method not in df["method"].values:
            continue
        grp = df[df["method"] == method].copy()
        results[method] = aggregate_method(grp, method)
        r = results[method]
        print(f"  {method}: passed={r['passed_all']:.1f} / {r['total_sequence_equivalents']:.0f}"
              f"  ({r['passed_all_rate']:.3f})")

    summary = {"categories": FAILURE_CATEGORIES, "methods": results}
    (out_dir / "data/failure_attribution_summary.json").write_text(json.dumps(summary, indent=2))
    print("  failure_attribution_summary.json saved")

    # Constraint overlap matrix
    overlap = compute_overlap_matrix(df)
    overlap.to_csv(out_dir / "data/constraint_overlap_matrix.csv", index=False)
    print("  constraint_overlap_matrix.csv saved")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="benchmarks/results/v3.2.0/benchmark_results.csv")
    parser.add_argument("--output-dir", default="reproducibility/benchmark_v0.5.1")
    args = parser.parse_args()
    run(Path(args.input), Path(args.output_dir))


if __name__ == "__main__":
    main()
