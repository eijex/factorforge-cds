# factorforge/benchmark/plot_paper5_figures.py
"""Paper 5 Publication Figures Generator.

Generates:
  Figure 2: The Biological Trilemma: Empirical Trade-offs & Sequential Greedy Distortion Trajectory
  Figure 3: Comparative Benchmark Across 5 Solver Arms (Feasibility, MFE, Runtime)
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = Path("c:/Work/eijex/eijex-workspace/_papers/manuscripts/_prism/paper5-factorforge-assembly/data")
FIGURES_DIR = Path("c:/Work/eijex/eijex-workspace/_papers/manuscripts/_prism/paper5-factorforge-assembly/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_CSV = DATA_DIR / "paper5_benchmark_results.csv"
TRAJECTORY_CSV = DATA_DIR / "greedy_distortion_trajectory.csv"
FIG2_PATH = FIGURES_DIR / "fig2_empirical_trilemma.png"
FIG3_PATH = FIGURES_DIR / "fig3_benchmark_comparison.png"


def load_results():
    rows = []
    with open(RESULTS_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            r["length_aa"] = int(r["length_aa"])
            r["cai_body"] = float(r["cai_body"])
            r["initiation_mfe_kcal"] = float(r["initiation_mfe_kcal"])
            r["bsai_sites"] = int(r["bsai_sites"])
            r["bsmbi_sites"] = int(r["bsmbi_sites"])
            r["runtime_ms"] = float(r["runtime_ms"])
            r["feasible"] = r["feasible"].lower() == "true"
            rows.append(r)
    return rows


def load_trajectory():
    rows = []
    with open(TRAJECTORY_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            r["step"] = int(r["step"])
            r["cai"] = float(r["cai"])
            r["mfe"] = float(r["mfe"])
            r["bsai"] = int(r["bsai"])
            r["bsmbi"] = int(r["bsmbi"])
            rows.append(r)
    return rows


def plot_figure_2(results, trajectory):
    """Plot Figure 2: The Biological Trilemma & Sequential Greedy Distortion Trajectory."""
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Panel A: Delta CAI vs Delta MFE Pareto Trade-off Plane
    dp21 = [r for r in results if r["solver_arm"] == "Arm_A_DP_v2_1"]
    dp20 = [r for r in results if r["solver_arm"] == "Arm_B_DP_v2_0"]
    greedy = [r for r in results if r["solver_arm"] == "Arm_C_Greedy"]
    codon = [r for r in results if r["solver_arm"] == "Arm_E_CoDOn"]

    ax1.scatter([r["cai_body"] for r in dp20], [r["initiation_mfe_kcal"] for r in dp20],
                color="#e74c3c", alpha=0.7, s=70, label="DP v2.0 (High CAI, Tight 5' MFE Bottleneck)")
    ax1.scatter([r["cai_body"] for r in dp21], [r["initiation_mfe_kcal"] for r in dp21],
                color="#2ecc71", alpha=0.9, s=90, edgecolors="#1b9e4b", linewidth=1.5,
                label="FactorForge DP v2.1 (Joint Exact Global Optimum)")
    ax1.scatter([r["cai_body"] for r in greedy], [r["initiation_mfe_kcal"] for r in greedy],
                color="#f39c12", alpha=0.6, s=60, marker="^", label="Sequential Greedy Patching")
    ax1.scatter([r["cai_body"] for r in codon], [r["initiation_mfe_kcal"] for r in codon],
                color="#9b59b6", alpha=0.6, s=60, marker="s", label="CoDOn (NSGA-II Genetic Algorithm)")

    ax1.axhline(y=-3.0, color="gray", linestyle="--", alpha=0.7, label="Threshold (MFE >= -3.0 kcal/mol)")
    ax1.set_xlabel("Body Codon Adaptation Index (CAI)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("5' Initiation Window MFE ($\Delta G$, kcal/mol)", fontsize=12, fontweight="bold")
    ax1.set_title("Panel A: Empirical Trilemma & Non-Dominated Frontier (N=36)", fontsize=13, fontweight="bold")
    ax1.legend(loc="lower left", fontsize=9, frameon=True)

    # Panel B: Sequential Greedy Distortion Trajectory
    steps = [0, 1, 2, 3, 4]
    step_labels = ["Step 0\n(Native)", "Step 1\n(CAI Max)", "Step 2\n(BsaI Fix)", "Step 3\n(BsmBI Fix)", "Step 4\n(GC Fix)"]
    
    mean_cai_per_step = [np.mean([t["cai"] for t in trajectory if t["step"] == s]) for s in steps]
    mean_mfe_per_step = [np.mean([t["mfe"] for t in trajectory if t["step"] == s]) for s in steps]

    color = "#2980b9"
    ax2.plot(steps, mean_cai_per_step, marker="o", linewidth=2.5, color=color, label="Mean CAI")
    ax2.set_xlabel("Sequential Heuristic Repair Pipeline", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Body CAI", color=color, fontsize=12, fontweight="bold")
    ax2.tick_params(axis="y", labelcolor=color)
    ax2.set_xticks(steps)
    ax2.set_xticklabels(step_labels, fontsize=10)

    ax2_twin = ax2.twinx()
    color_mfe = "#c0392b"
    ax2_twin.plot(steps, mean_mfe_per_step, marker="s", linewidth=2.5, linestyle="--", color=color_mfe, label="Mean 5' MFE")
    ax2_twin.set_ylabel("5' Initiation MFE ($\Delta G$, kcal/mol)", color=color_mfe, fontsize=12, fontweight="bold")
    ax2_twin.tick_params(axis="y", labelcolor=color_mfe)

    ax2.set_title("Panel B: Sequential Greedy Patching Distortion Trajectory", fontsize=13, fontweight="bold")
    fig.tight_layout()
    plt.savefig(FIG2_PATH, dpi=300)
    plt.close()
    print(f"Figure 2 saved to {FIG2_PATH}")


def plot_figure_3(results):
    """Plot Figure 3: Head-to-Head Benchmark Across 5 Solver Arms."""
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5), dpi=300)

    arm_keys = ["Arm_A_DP_v2_1", "Arm_B_DP_v2_0", "Arm_C_Greedy", "Arm_D_DNA_Chisel", "Arm_E_CoDOn"]
    arm_labels = ["DP v2.1\n(Exact)", "DP v2.0\n(Legacy)", "Sequential\nGreedy", "DNA\nChisel", "CoDOn\n(NSGA-II)"]
    palette = ["#2ecc71", "#e74c3c", "#f39c12", "#3498db", "#9b59b6"]

    # Panel A: 5' Initiation Window MFE Relaxation Boxplot
    mfe_data = [[r["initiation_mfe_kcal"] for r in results if r["solver_arm"] == k] for k in arm_keys]
    bp = ax1.boxplot(mfe_data, patch_artist=True, notch=True)
    ax1.set_xticks(range(1, len(arm_labels) + 1))
    ax1.set_xticklabels(arm_labels, fontsize=10)
    for patch, color in zip(bp["boxes"], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax1.axhline(y=-3.0, color="red", linestyle=":", alpha=0.8, label="Open Threshold (-3.0)")
    ax1.set_ylabel("5' Initiation Window MFE ($\Delta G$, kcal/mol)", fontsize=11, fontweight="bold")
    ax1.set_title("Panel A: 5' Initiation Structure Relaxation", fontsize=12, fontweight="bold")
    ax1.legend(loc="lower left", fontsize=9)

    # Panel B: Hard-Constraint Feasibility Rate (%)
    feas_rates = [
        sum(1 for r in results if r["solver_arm"] == k and r["feasible"]) / sum(1 for r in results if r["solver_arm"] == k) * 100.0
        for k in arm_keys
    ]
    bars = ax1_bars = ax2.bar(arm_labels, feas_rates, color=palette, alpha=0.85, edgecolor="black", linewidth=1.2)
    ax2.set_ylabel("Hard-Constraint Feasible Rate (%)", fontsize=11, fontweight="bold")
    ax2.set_ylim(0, 115)
    ax2.set_title("Panel B: Assembly Invariant Feasibility (BsaI/BsmBI=0)", fontsize=12, fontweight="bold")
    for bar in bars:
        yval = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2.0, yval + 2, f"{yval:.1f}%", ha="center", va="bottom", fontweight="bold")

    # Panel C: Runtime vs Sequence Length Scaling O(N)
    for k, label, color in zip(arm_keys, ["DP v2.1", "DP v2.0", "Greedy", "DNA Chisel", "CoDOn"], palette):
        arm_r = sorted([r for r in results if r["solver_arm"] == k], key=lambda x: x["length_aa"])
        lengths = [r["length_aa"] for r in arm_r]
        runtimes = [r["runtime_ms"] for r in arm_r]
        ax3.plot(lengths, runtimes, marker="o", label=label, color=color, linewidth=2, alpha=0.85)

    ax3.set_xlabel("Protein Length (Amino Acids)", fontsize=11, fontweight="bold")
    ax3.set_ylabel("Wall-Clock Runtime (ms)", fontsize=11, fontweight="bold")
    ax3.set_yscale("log")
    ax3.set_title("Panel C: Computational Runtime Scaling $O(N)$", fontsize=12, fontweight="bold")
    ax3.legend(loc="upper left", fontsize=9, frameon=True)

    fig.tight_layout()
    plt.savefig(FIG3_PATH, dpi=300)
    plt.close()
    print(f"Figure 3 saved to {FIG3_PATH}")


def main():
    print(f"Loading benchmark results from {RESULTS_CSV}...")
    results = load_results()
    trajectory = load_trajectory()
    plot_figure_2(results, trajectory)
    plot_figure_3(results)
    print("All Paper 5 figures generated successfully!")


if __name__ == "__main__":
    main()
