"""
Generate Figure 2, Figure 3, and Table 3 from benchmark_summary.json.
Does NOT require benchmark_results.csv.
Run from repo root:
  python reproducibility/benchmark_v0.5.1/scripts/figures/make_benchmark_figures.py
"""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Allow running from repo root without install
sys.path.insert(0, str(Path(__file__).parent.parent))
from common import METHOD_ORDER, DISPLAY_NAMES, BAR_COLORS

HEATMAP_METRICS = [
    ("pass_rate_multi_constraint", "Multi-constraint\npass rate"),
    ("mean_cai", "Mean CAI"),
    ("gc_in_range_rate", "GC in\ntarget range"),
]


def load_summary(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def make_figure2(summary: dict, out_dir: Path) -> None:
    methods_data = summary["methods"]
    pass_rates = [methods_data[m]["pass_rate_multi_constraint"] for m in METHOD_ORDER]
    colors = [BAR_COLORS[m] for m in METHOD_ORDER]
    labels = [DISPLAY_NAMES[m] for m in METHOD_ORDER]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(range(len(METHOD_ORDER)), pass_rates, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(METHOD_ORDER)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Multi-constraint pass rate", fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.set_title("Multi-constraint pass rate by design method\n(N. benthamiana, n=49,257)", fontsize=12, pad=12)
    ax.axhline(0.5, color="#9E9E9E", linewidth=0.5, linestyle="--")

    for bar, rate in zip(bars, pass_rates):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{rate:.3f}",
            ha="center", va="bottom", fontsize=8,
        )

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#BDBDBD", label="Baseline"),
        Patch(facecolor="#1976D2", label="FactorForge profile"),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc="upper left")
    fig.tight_layout()

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".svg"):
        fig.savefig(fig_dir / f"figure2_multiconstraint_pass_rate{suffix}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_figure3(summary: dict, out_dir: Path) -> None:
    methods_data = summary["methods"]
    metric_keys = [m[0] for m in HEATMAP_METRICS]
    metric_labels = [m[1] for m in HEATMAP_METRICS]
    labels = [DISPLAY_NAMES[m] for m in METHOD_ORDER]

    data = np.array([[methods_data[m][k] for k in metric_keys] for m in METHOD_ORDER])

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(data, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(HEATMAP_METRICS)))
    ax.set_xticklabels(metric_labels, fontsize=9)
    ax.set_yticks(range(len(METHOD_ORDER)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title("Benchmark metrics by design method", fontsize=12, pad=12)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for i in range(len(METHOD_ORDER)):
        for j in range(len(HEATMAP_METRICS)):
            val = data[i, j]
            color = "white" if val > 0.6 else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=8, color=color)

    fig.tight_layout()
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".svg"):
        fig.savefig(fig_dir / f"figure3_benchmark_tradeoff_heatmap{suffix}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_table3(summary: dict, out_dir: Path) -> None:
    methods_data = summary["methods"]
    headers = ["Method", "N", "Multi-constraint pass", "Mean CAI", "GC in range", "Biological pass", "Assembly pass"]
    rows = []
    for m in METHOD_ORDER:
        d = methods_data[m]
        rows.append([
            DISPLAY_NAMES[m].replace("\n", " "),
            str(d["n_ok"]),
            f"{d['pass_rate_multi_constraint']:.4f}",
            f"{d['mean_cai']:.4f}",
            f"{d['gc_in_range_rate']:.4f}",
            f"{d['pass_rate_biological']:.4f}",
            f"{d['pass_rate_assembly']:.4f}",
        ])

    sep = "| " + " | ".join(["---"] * len(headers)) + " |"
    lines = ["| " + " | ".join(headers) + " |", sep]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")

    tables_dir = out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    (tables_dir / "table3_benchmark_summary.md").write_text("\n".join(lines) + "\n")


def freeze_summary(summary_path: Path, out_dir: Path) -> None:
    sha256 = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    data_dir = out_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary_path, data_dir / "benchmark_summary.frozen.json")
    (data_dir / "benchmark_summary_sha256.json").write_text(
        json.dumps({"source": str(summary_path), "sha256": sha256}, indent=2)
    )


def main():
    parser = argparse.ArgumentParser(description="Generate manuscript figures from benchmark_summary.json")
    parser.add_argument("--input", default="benchmarks/results/v3.2.0/benchmark_summary.json")
    parser.add_argument("--output-dir", default="reproducibility/benchmark_v0.5.1")
    args = parser.parse_args()

    summary_path = Path(args.input)
    if not summary_path.exists():
        print(f"ERROR: {summary_path} not found. Run from repo root.", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.output_dir)
    summary = load_summary(summary_path)
    print(f"FactorForge {summary['factorforge_version']} | n={summary['dataset_n']} | seed={summary['random_seed']}")

    make_figure2(summary, out_dir)
    print("  Figure 2 saved")
    make_figure3(summary, out_dir)
    print("  Figure 3 saved")
    make_table3(summary, out_dir)
    print("  Table 3 saved")
    freeze_summary(summary_path, out_dir)
    print("  Summary frozen")
    print("Done.")


if __name__ == "__main__":
    main()
