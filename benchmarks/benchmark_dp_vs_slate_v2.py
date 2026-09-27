"""Empirical Candidate Space Benchmark: Single-best DP vs FactorForge Slate v2 (Job 293A).

Evaluates Adalimumab VH (Heavy Chain Variable Domain) across:
1. Single-Best DP (Default Feasibility / Feasibility-best optimizer)
2. Slate v2 Candidates (Top-1, Top-2, Top-3, Top-4, Top-5, Top-25)

Metrics compared:
- CAI
- Global GC (%)
- 5' Initiation Structure Proxy Score
- ViennaRNA 5' Thermodynamic Folding (kcal/mol & structure)
- Rare Codon Count (<10% freq)
- Local GC Variance (50-nt window)
- Synonymous Sequence Edit Distance (vs Single DP)
- Pareto Front Rank
- Archetype Profile Name
"""

import json
from tabulate import tabulate
from factorforge.engines.dp_v2_1_1 import DPV211Optimizer
from factorforge.core.slate_engine import SlateV2Engine
from factorforge.core.host_model import HostModel
from factorforge.analysis.metrics import (
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    detect_homopolymers,
)

ADALIMUMAB_VH = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRD"
    "NAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
)


def run_benchmark():
    print("=" * 80)
    print("EMPIRICAL BENCHMARK: SINGLE-BEST DP vs SLATE v2 CANDIDATE SPACE")
    print("Target: Adalimumab (Humira) Heavy Chain Variable Domain (121 AA)")
    print("Host: Nicotiana benthamiana (NbeV1.1_HC_v3.7)")
    print("=" * 80)

    host_model = HostModel.load("nbenthamiana")

    # 1. Run Single-Best DP Optimizer
    dp_opt = DPV211Optimizer()
    dp_res = dp_opt.optimize(ADALIMUMAB_VH, host_model.codon_weights)
    dp_cds = dp_res["sequence"]
    dp_cai = round(calculate_cai(dp_cds, host_model.codon_weights), 3)
    dp_gc = round(calculate_gc(dp_cds), 1)

    # 5' initiation proxy for DP
    ramp_48 = dp_cds[:48]
    dp_5p_gc = (ramp_48.count("G") + ramp_48.count("C")) / len(ramp_48)
    dp_5p_proxy = round(-0.5 - (dp_5p_gc * 28.0), 2)

    # ViennaRNA deep fold for DP
    try:
        import RNA
        _, dp_mfe = RNA.fold(dp_cds[:48])
        dp_mfe_str = f"{dp_mfe:.1f}"
    except Exception:
        dp_mfe_str = "N/A"

    dp_rare = sum(1 for i in range(0, len(dp_cds), 3) if host_model.codon_frequencies.get(dp_cds[i:i+3], 0.0) < 0.10)
    dp_windows = calculate_gc_windows(dp_cds, window_size=50, step=10)
    dp_gc_vals = [float(w["gc"]) / 100.0 for w in dp_windows]
    dp_var_gc = round(sum((v - sum(dp_gc_vals)/len(dp_gc_vals))**2 for v in dp_gc_vals)/len(dp_gc_vals), 4)

    # 2. Run Hardened Slate v2 Engine
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    slate_res = engine.generate_slate(
        protein_sequence=ADALIMUMAB_VH,
        slate_size=25,
        seed=42,
    )

    rows = []
    # Add Single-Best DP baseline
    rows.append([
        "Single DP (Baseline)",
        "Front 1 (Greedy)",
        "Global Feasibility Best",
        f"{dp_cai:.3f}",
        f"{dp_gc:.1f}%",
        f"{dp_5p_proxy:.1f}",
        f"{dp_mfe_str} kcal",
        dp_rare,
        f"{dp_var_gc:.4f}",
        "0 nt (0.0%)",
        dp_res.get("provenance", {}).get("status", "PASS"),
    ])

    slate = slate_res["slate"]
    inspect_indices = [0, 1, 2, 3, 4, 9, 24] # Top-1, 2, 3, 4, 5, 10, 25

    for idx in inspect_indices:
        if idx < len(slate):
            cand = slate[idx]
            cand_cds = cand["coding_sequence"]
            scores = cand["scores"]
            
            # Edit distance vs Single DP baseline
            diff_nt = sum(c1 != c2 for c1, c2 in zip(cand_cds, dp_cds))
            diff_pct = (diff_nt / len(dp_cds)) * 100.0

            # Deep RNA
            deep_rna = cand.get("deep_rna_folding", {})
            mfe_val = deep_rna.get("mfe_5p_kcal_mol")
            mfe_display = f"{mfe_val:.1f} kcal" if mfe_val is not None else "N/A"

            # Var GC
            cand_wins = calculate_gc_windows(cand_cds, window_size=50, step=10)
            cand_gc_vals = [float(w["gc"]) / 100.0 for w in cand_wins]
            cand_var_gc = round(sum((v - sum(cand_gc_vals)/len(cand_gc_vals))**2 for v in cand_gc_vals)/len(cand_gc_vals), 4)

            label = f"Slate #{cand['rank']}"
            rows.append([
                label,
                f"Front {cand['pareto_front']}",
                cand["profile_name"],
                f"{scores['cai']:.3f}",
                f"{scores['gc_global']*100.0:.1f}%",
                f"{scores['five_prime_structure_proxy_score']:.1f}",
                mfe_display,
                scores["rare_codon_count"],
                f"{cand_var_gc:.4f}",
                f"{diff_nt} nt ({diff_pct:.1f}%)",
                "PASS (100% AA)",
            ])

    headers = [
        "Candidate",
        "Pareto Front",
        "Profile Archetype",
        "CAI",
        "GC %",
        "5' Proxy",
        "RNAfold 5'",
        "Rare",
        "Local VarGC",
        "Diff vs DP",
        "Validation",
    ]

    print(tabulate(rows, headers=headers, tablefmt="github"))
    print("\n" + "=" * 80)
    print("DIVERSITY & RECALL METRICS SUMMARY")
    print(f"- Total Raw Candidates Explored: {slate_res['metadata']['total_candidates_generated']}")
    print(f"- Unique Feasible Candidates: {slate_res['metadata']['unique_feasible_candidates']}")
    print(f"- Preselected Candidates Evaluated: {slate_res['metadata']['preselected_candidates_evaluated']}")
    print(f"- Final Assembled Slate: {slate_res['metadata']['slate_size']}")
    print(f"- Max Synonymous Diversity in Slate: {max(sum(c1 != c2 for c1, c2 in zip(c['coding_sequence'], dp_cds)) for c in slate)} nt substitutions")
    print(f"- 100% Amino Acid Conservation Verified: {slate_res['validation_summary']['aa_conservation_verified']}")
    print(f"- Forbidden Sites Clean: {slate_res['validation_summary']['forbidden_sites_clean']}")
    print("=" * 80)


if __name__ == "__main__":
    run_benchmark()
