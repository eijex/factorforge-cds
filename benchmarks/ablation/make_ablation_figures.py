# benchmarks/ablation/make_ablation_figures.py
"""Generate ablation figures from ablation_summary.json.

Usage:
    python -m benchmarks.ablation.make_ablation_figures
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SUMMARY_PATH = ROOT / "benchmarks" / "results" / "v3.2.0" / "ablation" / "ablation_summary.json"
OUT_DIR = ROOT / "benchmarks" / "results" / "v3.2.0" / "ablation" / "figures"

LAYER_ORDER = ["L0", "L1", "L2", "L3", "L4", "L5"]
LAYER_LABELS = {
    "L0": "L0\nrandom_syn.",
    "L1": "L1\ngreedy_cai",
    "L2": "L2\ncai+gc",
    "L3": "L3\ncai+iis",
    "L4": "L4\ncai+gc+iis",
    "L5": "L5\nassembly_friendly",
}
BAR_COLORS = {
    "L0": "#B0BEC5",
    "L1": "#64B5F6",
    "L2": "#FFB300",
    "L3": "#FB8C00",
    "L4": "#E53935",
    "L5": "#43A047",
}
HEATMAP_METRICS = [
    ("multi_constraint_pass_rate", "Multi-constraint\npass rate"),
    ("mean_cai", "Mean CAI"),
    ("gc_in_range_rate", "GC in-range rate"),
    ("assembly_pass_rate", "Assembly pass rate"),
]


def make_ablation_pass_rate_figure(summary: dict, out_dir: Path) -> None:
    layers = summary["layers"]
    present = [l for l in LAYER_ORDER if l in layers]
    rates = [layers[l]["multi_constraint_pass_rate"] for l in present]
    labels = [LAYER_LABELS[l] for l in present]
    colors = [BAR_COLORS[l] for l in present]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(present))
    bars = ax.bar(x, rates, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Multi-constraint pass rate", fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
    ax.set_title("L0–L5 Ablation: Multi-constraint Pass Rate", fontsize=11, fontweight="bold")
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{rate:.1%}", ha="center", va="bottom", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "figure_ablation_pass_rate.png", dpi=150)
    fig.savefig(out_dir / "figure_ablation_pass_rate.svg")
    plt.close(fig)


def make_ablation_heatmap_figure(summary: dict, out_dir: Path) -> None:
    layers = summary["layers"]
    present = [l for l in LAYER_ORDER if l in layers]
    metric_keys = [k for k, _ in HEATMAP_METRICS]
    metric_labels = [lbl for _, lbl in HEATMAP_METRICS]

    data = np.array([[layers[l].get(m, 0.0) for m in metric_keys] for l in present])
    row_labels = [LAYER_LABELS[l] for l in present]

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(data, aspect="auto", cmap="YlGn", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(metric_keys)))
    ax.set_xticklabels(metric_labels, fontsize=9)
    ax.set_yticks(np.arange(len(present)))
    ax.set_yticklabels(row_labels, fontsize=9)
    for i in range(len(present)):
        for j in range(len(metric_keys)):
            val = data[i, j]
            color = "white" if val > 0.65 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color=color)
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    ax.set_title("L0–L5 Ablation: Constraint Trade-off Heatmap", fontsize=11, fontweight="bold")
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "figure_ablation_tradeoff_heatmap.png", dpi=150)
    fig.savefig(out_dir / "figure_ablation_tradeoff_heatmap.svg")
    plt.close(fig)


def main() -> None:
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    make_ablation_pass_rate_figure(summary, OUT_DIR)
    print(f"[ablation-figures] pass rate figure → {OUT_DIR}/figure_ablation_pass_rate.png")
    make_ablation_heatmap_figure(summary, OUT_DIR)
    print(f"[ablation-figures] heatmap figure → {OUT_DIR}/figure_ablation_tradeoff_heatmap.png")


if __name__ == "__main__":
    main()
