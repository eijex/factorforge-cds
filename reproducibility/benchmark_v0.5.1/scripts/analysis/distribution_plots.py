"""
Distribution plots for CAI and GC% from benchmark_results.csv.
random_synonymous: all 3 replicates plotted as separate points.
Run from repo root:
  python reproducibility/benchmark_v0.5.1/scripts/analysis/distribution_plots.py
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from common import METHOD_ORDER, DISPLAY_NAMES, BAR_COLORS


def plot_distribution(df: pd.DataFrame, column: str, xlabel: str,
                      title: str, filename: str, out_dir: Path) -> None:
    """Violin plot per method. random_synonymous: all 3 replicates included."""
    fig, axes = plt.subplots(1, 1, figsize=(12, 5))
    ax = axes

    positions = list(range(len(METHOD_ORDER)))
    data_per_method = []
    labels = []
    for m in METHOD_ORDER:
        grp = df[df["method"] == m][column].dropna().astype(float)
        data_per_method.append(grp.values)
        labels.append(DISPLAY_NAMES[m])

    vp = ax.violinplot(data_per_method, positions=positions, showmedians=True,
                       widths=0.7)
    for i, (body, m) in enumerate(zip(vp["bodies"], METHOD_ORDER)):
        body.set_facecolor(BAR_COLORS[m])
        body.set_alpha(0.5)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel(xlabel, fontsize=11)
    ax.set_title(title, fontsize=12)

    fig.tight_layout()
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {filename}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="benchmarks/results/v3.2.0/benchmark_results.csv")
    parser.add_argument("--output-dir", default="reproducibility/benchmark_v0.5.1")
    args = parser.parse_args()

    print(f"Loading {args.input} (CAI and gc_percent columns only) ...")
    df = pd.read_csv(args.input, usecols=["method", "sequence_id", "replicate", "cai", "gc_percent"])
    out_dir = Path(args.output_dir)

    plot_distribution(
        df, "cai", "CAI",
        "CAI distribution by design method (random_synonymous: all 3 replicates)",
        "figure_cai_distribution.png", out_dir,
    )
    plot_distribution(
        df, "gc_percent", "GC content (%)",
        "GC% distribution by design method",
        "figure_gc_distribution.png", out_dir,
    )
    print("Done.")


if __name__ == "__main__":
    main()
