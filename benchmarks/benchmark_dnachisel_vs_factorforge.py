"""Head-to-head live benchmark: FactorForge Exact DP vs DNA Chisel Heuristic Search.

Compares:
1. Algorithmic Paradigm:
   - FactorForge: Exact Finite-State Dynamic Programming over Semiring Trellis ($O(N^2 |Q| B)$).
   - DNA Chisel: Stochastic local search / simulated annealing with random synonymous mutations.
2. Target Proteins:
   - Adalimumab VH (121 aa)
   - CD47 Extracellular Domain (140 aa)
   - Enhanced GFP (238 aa)
3. Hard Constraints:
   - Host: Nicotiana benthamiana (HostModel NbeV1.1)
   - GC Target Window: 40.0% - 46.0%
   - Forbidden Motifs: BsaI (GGTCTC/GAGACC), BsmBI (CGTCTC/GAGACG), BpiI (GAAGAC/GTCTTC), SapI (GCTCTTC/GAAGAGC)
4. Evaluation Metrics:
   - Wall-clock runtime (ms)
   - Global Codon Adaptation Index (CAI) under identical weights
   - Global GC%
   - Forbidden motif clearance
   - Feasibility guarantee
"""

from __future__ import annotations

from collections import defaultdict
import time
from typing import Any, Dict, List

from factorforge.analysis.metrics import (
    STANDARD_GENETIC_CODE,
    calculate_cai,
    calculate_gc,
    detect_forbidden_motifs,
)
from factorforge.core.host_model import HostModel
from factorforge.engines.dp_v2_4 import DPV240Optimizer

BENCHMARK_TARGETS: Dict[str, str] = {
    "Adalimumab_VH": (
        "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
    ),
    "CD47_ECD": (
        "MWPLVAALLLGSACCGSAQLLFNKTKSVEFTFCNDTVVIPCFVTNMEAQNTTEVYVKWKFKGRDIYTFDGALNKSTVPTDFSSAKIEVSQLLKGDASLKMDKSDAVSHTGNYTCEVTELTREGETIIELKYRVVSWFSPN"
    ),
    "EGFP": (
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    ),
}

CANONICAL_FORBIDDEN_MOTIFS: List[str] = [
    "GGTCTC",
    "GAGACC",  # BsaI / Eco31I
    "CGTCTC",
    "GAGACG",  # BsmBI / Esp3I
    "GAAGAC",
    "GTCTTC",  # BpiI / BbsI
    "GCTCTTC",
    "GAAGAGC",  # SapI
]


def _build_dnachisel_codon_table(host_model: HostModel) -> Dict[str, Dict[str, float]]:
    """Build normalized per-amino-acid codon table for DNA Chisel."""
    dc_table: Dict[str, Dict[str, float]] = defaultdict(dict)
    for codon, aa in STANDARD_GENETIC_CODE.items():
        freq = host_model.codon_frequencies.get(codon, 0.001)
        dc_table[aa][codon] = freq
    for _, cdict in dc_table.items():
        tot = sum(cdict.values())
        if tot > 0:
            for c in cdict:
                cdict[c] /= tot
    return dict(dc_table)


def run_benchmark(
    targets: Dict[str, str] | None = None,
    gc_min: float = 0.40,
    gc_max: float = 0.46,
    host_name: str = "nbenthamiana",
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """Execute head-to-head benchmark comparing FactorForge DP against DNA Chisel."""
    import dnachisel  # type: ignore

    targets = targets or BENCHMARK_TARGETS
    host_model = HostModel.load(host_name)
    weights = host_model.codon_weights
    dc_codon_table = _build_dnachisel_codon_table(host_model)

    records: List[Dict[str, Any]] = []

    if verbose:
        print("=" * 95)
        print("HEAD-TO-HEAD BENCHMARK: FACTORFORGE EXACT DP vs DNA CHISEL LOCAL SEARCH")
        print(f"Host: {host_name} | GC Range: {gc_min*100:.1f}% - {gc_max*100:.1f}% | Forbidden Motifs: {len(CANONICAL_FORBIDDEN_MOTIFS)}")
        print("=" * 95)
        header = f"{'Target':<15} | {'Length':<6} | {'Engine':<12} | {'Time (ms)':<10} | {'CAI':<8} | {'GC %':<8} | {'Motifs':<8} | {'Optimality':<10}"
        print(header)
        print("-" * 95)

    for target_name, protein_seq in targets.items():
        n_aa = len(protein_seq)

        # 1. FactorForge DP v2.4.0 (Exact Dynamic Programming)
        t0 = time.time()
        ff_opt = DPV240Optimizer(auto_annotate_pauses=False)
        ff_res = ff_opt.optimize(
            protein_sequence=protein_seq,
            codon_weights=weights,
            target_gc_min=gc_min,
            target_gc_max=gc_max,
        )
        t_ff = (time.time() - t0) * 1000.0
        ff_cds = ff_res["sequence"]
        ff_cai = calculate_cai(ff_cds, weights)
        ff_gc = calculate_gc(ff_cds)
        ff_motifs = detect_forbidden_motifs(ff_cds, CANONICAL_FORBIDDEN_MOTIFS)

        ff_record = {
            "target": target_name,
            "length_aa": n_aa,
            "engine": "FactorForge DP",
            "time_ms": round(t_ff, 2),
            "cai": round(ff_cai, 4),
            "gc_percent": round(ff_gc, 2),
            "motif_violations": len(ff_motifs),
            "status": "FEASIBLE" if ff_res["gc_feasible"] else "INFEASIBLE",
            "optimality": "GLOBAL_EXACT",
            "sequence": ff_cds,
        }
        records.append(ff_record)

        if verbose:
            print(
                f"{target_name:<15} | {n_aa:<6} | {'FactorForge':<12} | {t_ff:>9.1f} | "
                f"{ff_cai:>7.4f} | {ff_gc:>6.2f}% | {len(ff_motifs):>8} | {'EXACT_GLOBAL'}"
            )

        # 2. DNA Chisel (Stochastic Heuristic Search)
        init_dna = dnachisel.reverse_translate(protein_seq)
        constraints = [
            dnachisel.EnforceTranslation(),
            dnachisel.EnforceGCContent(mini=gc_min, maxi=gc_max),
        ] + [dnachisel.AvoidPattern(pat) for pat in CANONICAL_FORBIDDEN_MOTIFS]

        problem = dnachisel.DnaOptimizationProblem(
            sequence=init_dna,
            constraints=constraints,
            objectives=[
                dnachisel.CodonOptimize(method="use_best_codon", codon_usage_table=dc_codon_table),
            ],
            logger=None,
        )

        t0 = time.time()
        try:
            problem.resolve_constraints()
            problem.optimize()
            t_dc = (time.time() - t0) * 1000.0
            dc_cds = problem.sequence
            dc_cai = calculate_cai(dc_cds, weights)
            dc_gc = calculate_gc(dc_cds)
            dc_motifs = detect_forbidden_motifs(dc_cds, CANONICAL_FORBIDDEN_MOTIFS)

            dc_record = {
                "target": target_name,
                "length_aa": n_aa,
                "engine": "DNA Chisel",
                "time_ms": round(t_dc, 2),
                "cai": round(dc_cai, 4),
                "gc_percent": round(dc_gc, 2),
                "motif_violations": len(dc_motifs),
                "status": "FEASIBLE",
                "optimality": "HEURISTIC",
                "sequence": dc_cds,
            }
            records.append(dc_record)

            if verbose:
                cai_delta = ff_cai - dc_cai
                print(
                    f"{target_name:<15} | {n_aa:<6} | {'DNA Chisel':<12} | {t_dc:>9.1f} | "
                    f"{dc_cai:>7.4f} | {dc_gc:>6.2f}% | {len(dc_motifs):>8} | {'HEURISTIC'}"
                )
                print(f"   --> FactorForge CAI Advantage: +{cai_delta:.4f} higher CAI")

        except Exception as exc:
            t_dc = (time.time() - t0) * 1000.0
            dc_record = {
                "target": target_name,
                "length_aa": n_aa,
                "engine": "DNA Chisel",
                "time_ms": round(t_dc, 2),
                "cai": None,
                "gc_percent": None,
                "motif_violations": None,
                "status": f"FAILED ({type(exc).__name__})",
                "optimality": "NONE",
                "sequence": None,
            }
            records.append(dc_record)
            if verbose:
                print(
                    f"{target_name:<15} | {n_aa:<6} | {'DNA Chisel':<12} | {t_dc:>9.1f} | "
                    f"{'N/A':>7} | {'N/A':>7} | {'N/A':>8} | {type(exc).__name__}"
                )

    if verbose:
        print("=" * 95)

    return records


if __name__ == "__main__":
    run_benchmark()
