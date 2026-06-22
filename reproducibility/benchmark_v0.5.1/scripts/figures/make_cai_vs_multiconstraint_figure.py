"""
Scatter plot of mean CAI vs. corrected multi-constraint pass rate per method,
from benchmark_summary.frozen.json (scoring_contract v1.1).

This script is a reconstruction, not a recovery: no generation script for
this figure ("CAI-only baseline versus constraint-aware design") was found
anywhere in this repository or its Git history. The data points below were
verified to match the figure's existing seven per-method values (mean_cai,
pass_rate_multi_constraint) before this script was written.

Run from repo root:
  python reproducibility/benchmark_v0.5.1/scripts/figures/make_cai_vs_multiconstraint_figure.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Order and display labels match the published figure's legend, not common.py's
# bar-chart DISPLAY_NAMES (this figure predates that convention).
METHOD_DISPLAY = {
    "native_reference": "Native",
    "random_synonymous": "Random",
    "factorforge_high_cai": "High-CAI",
    "factorforge_balanced": "Balanced",
    "factorforge_gc_target": "GC",
    "greedy_cai": "Greedy",
    "factorforge_assembly_friendly": "AF",
}

# Points whose markers sit close together in the published figure carry an
# explicit text-box annotation. The published figure used the "ff_" shorthand
# for FactorForge methods rather than the internal method id.
ANNOTATED = {"factorforge_assembly_friendly", "factorforge_gc_target", "greedy_cai"}
ANNOTATION_LABEL = {
    "factorforge_assembly_friendly": "ff_assembly_friendly",
    "factorforge_gc_target": "ff_gc_target",
    "greedy_cai": "greedy_cai",
}

COLORS = {
    "native_reference": "#616161",
    "random_synonymous": "#8D6E1A",
    "factorforge_high_cai": "#D81B60",
    "factorforge_balanced": "#5E35B1",
    "factorforge_gc_target": "#E64A19",
    "greedy_cai": "#558B2F",
    "factorforge_assembly_friendly": "#00897B",
}


def main():
    summary_path = Path("reproducibility/benchmark_v0.5.1/data/benchmark_summary.frozen.json")
    if not summary_path.exists():
        print(f"ERROR: {summary_path} not found.", file=sys.stderr)
        sys.exit(1)

    summary = json.load(open(summary_path))
    methods_data = summary["methods"]

    fig, ax = plt.subplots(figsize=(9, 6))

    for method, label in METHOD_DISPLAY.items():
        d = methods_data[method]
        x = d["mean_cai"]
        y = d["pass_rate_multi_constraint"] * 100
        ax.scatter(x, y, s=120, color=COLORS[method], label=label, zorder=3)
        if method in ANNOTATED:
            ax.annotate(
                ANNOTATION_LABEL[method],
                xy=(x, y),
                xytext=(x + 0.04, y + 4),
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.5"),
                arrowprops=dict(arrowstyle="-", color="0.4", lw=0.8),
            )

    ax.set_xlabel("Mean general CAI", fontsize=11)
    ax.set_ylabel("Corrected multi-constraint pass rate (%)", fontsize=11)
    ax.set_title("CAI-only baseline versus constraint-aware design", fontsize=13)
    ax.set_xlim(0.70, 1.03)
    ax.set_ylim(-3, 70)
    ax.grid(True, alpha=0.4)
    ax.legend(title="Methods", loc="upper left")

    fig.tight_layout()
    out = Path("reproducibility/benchmark_v0.5.1/figures/figure_cai_vs_multiconstraint.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

    for method, label in METHOD_DISPLAY.items():
        d = methods_data[method]
        print(f"{label:10s} ({method:30s}) mean_cai={d['mean_cai']:.4f} "
              f"pass_rate={100*d['pass_rate_multi_constraint']:.2f}%")


if __name__ == "__main__":
    main()
