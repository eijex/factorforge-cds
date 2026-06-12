"""
Stacked bar chart of failure categories per method from failure_attribution_summary.json.
Run from repo root:
  python reproducibility/benchmark_v0.5.1/scripts/figures/make_attribution_figure.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from common import METHOD_ORDER, DISPLAY_NAMES, FAILURE_CATEGORIES

CATEGORY_COLORS = {
    "passed_all": "#43A047",
    "failed_gc_only": "#FFB300",
    "failed_iis_only": "#FB8C00",
    "failed_gc_and_iis": "#E53935",
    "failed_bio": "#8E24AA",
    "failed_other": "#B0BEC5",
}

CATEGORY_LABELS = {
    "passed_all": "Passed all",
    "failed_gc_only": "Failed: GC only",
    "failed_iis_only": "Failed: Type IIS only",
    "failed_gc_and_iis": "Failed: GC + Type IIS",
    "failed_bio": "Failed: biological",
    "failed_other": "Failed: other",
}


def main():
    summary_path = Path("reproducibility/benchmark_v0.5.1/data/failure_attribution_summary.json")
    if not summary_path.exists():
        print(f"ERROR: {summary_path} not found. Run failure_attribution.py first.", file=sys.stderr)
        sys.exit(1)

    summary = json.load(open(summary_path))
    methods_data = summary["methods"]

    present = [m for m in METHOD_ORDER if m in methods_data]
    labels = [DISPLAY_NAMES[m] for m in present]
    n = len(present)

    fig, ax = plt.subplots(figsize=(11, 5))
    bottoms = np.zeros(n)

    for cat in FAILURE_CATEGORIES:
        rates = np.array([methods_data[m].get(f"{cat}_rate", 0.0) for m in present])
        ax.bar(range(n), rates, bottom=bottoms,
               color=CATEGORY_COLORS[cat], label=CATEGORY_LABELS[cat], edgecolor="white", linewidth=0.3)
        bottoms += rates

    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Proportion of sequences", fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.set_title("Failure attribution by design method\n(N. benthamiana, n=49,257)", fontsize=12, pad=12)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)

    fig.tight_layout()
    out = Path("reproducibility/benchmark_v0.5.1/figures/figure_failure_attribution.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
