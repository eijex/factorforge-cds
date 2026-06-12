"""
Sequence-level paired bootstrap CI for multi_constraint_pass rate.
random_synonymous: mean across replicates per sequence before resampling.
Run from repo root:
  python reproducibility/benchmark_v0.5.1/scripts/analysis/bootstrap_ci.py
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from common import METHOD_ORDER, REPLICATED_METHOD

DIFFERENCE_PAIRS = [
    ("factorforge_assembly_friendly", "greedy_cai"),
    ("factorforge_assembly_friendly", "random_synonymous"),
    ("factorforge_assembly_friendly", "native_reference"),
]


def build_sequence_pass_matrix(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """
    Returns {method: array of shape (n_sequences,)} with per-sequence pass rates.
    For random_synonymous: mean across replicates. For others: single value.
    All arrays are aligned to the same sequence order.
    """
    sequences = sorted(df["sequence_id"].unique())
    seq_index = {s: i for i, s in enumerate(sequences)}
    n = len(sequences)

    result = {}
    for method in METHOD_ORDER:
        grp = df[df["method"] == method]
        if len(grp) == 0:
            continue
        arr = np.zeros(n)
        if method == REPLICATED_METHOD and "replicate" in grp.columns and grp["replicate"].nunique() > 1:
            means = grp.groupby("sequence_id")["multi_constraint_pass"].mean()
            for seq, val in means.items():
                arr[seq_index[seq]] = float(val)
        else:
            idxs = grp["sequence_id"].map(seq_index).values
            arr[idxs] = grp["multi_constraint_pass"].astype(float).values
        result[method] = arr

    return result


def paired_bootstrap(
    seq_pass: dict[str, np.ndarray],
    n_resamples: int = 1000,
    seed: int = 320,
) -> dict:
    """Paired bootstrap over sequence_ids. Returns CI dict."""
    rng = np.random.default_rng(seed)
    n_seq = len(next(iter(seq_pass.values())))
    methods = list(seq_pass.keys())

    boot_rates = {m: np.zeros(n_resamples) for m in methods}
    for i in range(n_resamples):
        idx = rng.integers(0, n_seq, size=n_seq)
        for m, arr in seq_pass.items():
            boot_rates[m][i] = arr[idx].mean()

    per_method = {}
    for m, arr in seq_pass.items():
        rates = boot_rates[m]
        per_method[m] = {
            "observed_rate": float(arr.mean()),
            "ci_lower_2_5": float(np.percentile(rates, 2.5)),
            "ci_upper_97_5": float(np.percentile(rates, 97.5)),
        }

    difference_cis = {}
    af_rates = boot_rates.get("factorforge_assembly_friendly")
    if af_rates is not None:
        for _, ref in DIFFERENCE_PAIRS:
            if ref not in boot_rates:
                continue
            key = f"factorforge_assembly_friendly - {ref}"
            diff = af_rates - boot_rates[ref]
            obs_diff = (seq_pass["factorforge_assembly_friendly"].mean()
                        - seq_pass[ref].mean())
            difference_cis[key] = {
                "observed_difference": float(obs_diff),
                "ci_lower_2_5": float(np.percentile(diff, 2.5)),
                "ci_upper_97_5": float(np.percentile(diff, 97.5)),
            }

    return {
        "bootstrap_unit": "sequence_id",
        "paired_across_methods": True,
        "random_synonymous_aggregation": "mean_across_replicates_per_sequence",
        "n_resamples": n_resamples,
        "random_seed": seed,
        "metrics": ["multi_constraint_pass_rate"],
        "per_method": per_method,
        "difference_cis": difference_cis,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="benchmarks/results/v3.2.0/benchmark_results.csv")
    parser.add_argument("--output-dir", default="reproducibility/benchmark_v0.5.1")
    parser.add_argument("--n-resamples", type=int, default=1000)
    args = parser.parse_args()

    print(f"Loading {args.input} ...")
    df = pd.read_csv(args.input, usecols=[
        "method", "sequence_id", "replicate", "multi_constraint_pass"
    ])
    df["multi_constraint_pass"] = df["multi_constraint_pass"].astype(bool)

    print("Building sequence pass matrix ...")
    seq_pass = build_sequence_pass_matrix(df)
    print(f"  {len(seq_pass)} methods, {len(next(iter(seq_pass.values())))} sequences")

    print(f"Running bootstrap (n={args.n_resamples}) ...")
    result = paired_bootstrap(seq_pass, n_resamples=args.n_resamples, seed=320)

    out_path = Path(args.output_dir) / "data/bootstrap_ci.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2))
    print(f"  Saved to {out_path}")

    print("\nPer-method 95% CI:")
    for m, ci in result["per_method"].items():
        print(f"  {m}: {ci['observed_rate']:.4f} [{ci['ci_lower_2_5']:.4f}, {ci['ci_upper_97_5']:.4f}]")


if __name__ == "__main__":
    main()
