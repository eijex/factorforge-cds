"""Experimental 5DP Prototype: WSSV VP28-HBcAg VLP Assembly & Plant Glycan Layer.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Implements a research-preview 5DP scorer/generator for plant-expressed VLP
vaccines (WSSV VP28 ectodomain fused to HBcAg nanoscaffold in N. benthamiana
and Wolffia globosa), emitting wssv_yield_prediction_candidates.jsonl schemas.
"""

from __future__ import annotations

import json
import random
from collections.abc import Sequence
from typing import Any

from factorforge.analysis.metrics import calculate_cai, calculate_gc
from factorforge.engines.dp_v2_4 import DPV240Optimizer
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence

PR1A_SIGNAL_PEPTIDE = "MGFVLFSQLPSFLLVSTLLLFLVISHSCRA"
KDEL_TAG = "KDEL"
HIS6_TAG = "HHHHHH"
GGGGS_LINKER = "GGGGS"


def detect_n_glycosylation_motifs(protein: str) -> list[dict[str, Any]]:
    """Detect canonical N-X-S/T (X != P) N-glycosylation sequons in a protein."""
    seq = "".join(protein.upper().split()).rstrip("*")
    motifs: list[dict[str, Any]] = []
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in ("S", "T"):
            motifs.append(
                {
                    "position_1based": i + 1,
                    "sequon": seq[i : i + 3],
                    "risk": "PLANT_N_GLYCAN_SITE",
                }
            )
    return motifs


def assemble_wssv_vlp_construct(
    hbcag_seq: str,
    vp28_seq: str,
    topology: str = "c_terminal",
    include_pr1a: bool = True,
    include_kdel: bool = True,
    include_his6: bool = True,
    linker: str = GGGGS_LINKER,
) -> dict[str, Any]:
    """Assemble a WSSV VP28-HBcAg VLP fusion protein and compute boundary pause indices."""
    hbc = "".join(hbcag_seq.upper().split()).rstrip("*")
    vp = "".join(vp28_seq.upper().split()).rstrip("*")
    if not hbc or not vp:
        raise ValueError("Both hbcag_seq and vp28_seq must be non-empty.")

    prefix = PR1A_SIGNAL_PEPTIDE if include_pr1a else ""
    suffix = (HIS6_TAG if include_his6 else "") + (KDEL_TAG if include_kdel else "")

    topo = topology.strip().lower()
    if topo == "c_terminal":
        core = hbc + linker + vp
        boundary_aa = len(prefix) + len(hbc)
        steric_note = "C-terminal fusion; flexible linker separates HBcAg core from VP28."
    elif topo == "n_terminal":
        core = vp + linker + hbc
        boundary_aa = len(prefix) + len(vp)
        steric_note = "N-terminal fusion; VP28 precedes HBcAg assembly domain."
    elif topo == "mir_loop":
        mid = len(hbc) // 2
        core = hbc[:mid] + linker + vp + linker + hbc[mid:]
        boundary_aa = len(prefix) + mid
        steric_note = "MIR-loop insertion; surface-exposed VLP display (steric check required)."
    else:
        raise ValueError("topology must be 'c_terminal', 'n_terminal', or 'mir_loop'.")

    full_protein = prefix + core + suffix
    # 5DP boundary + exit-tunnel offset (+35 aa) pause window
    pause_center = min(len(full_protein), boundary_aa + 35)
    pause_positions = [p for p in (pause_center - 1, pause_center, pause_center + 1) if 1 <= p <= len(full_protein)]

    return {
        "topology": topo,
        "protein_sequence": full_protein,
        "length_aa": len(full_protein),
        "fusion_boundary_aa": boundary_aa,
        "pause_positions": pause_positions,
        "n_glycan_motifs": detect_n_glycosylation_motifs(vp),
        "has_pr1a": include_pr1a,
        "has_kdel": include_kdel,
        "steric_note": steric_note,
    }


class WSSV5DPOptimizer:
    """Experimental 5DP candidate generator for WSSV VP28-HBcAg VLP constructs."""

    def __init__(self, forbidden_motifs: list[str] | None = None) -> None:
        self.forbidden_motifs = forbidden_motifs or ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]

    def generate_candidates(
        self,
        hbcag_seq: str,
        vp28_seq: str,
        codon_weights: dict[str, float],
        host: str = "nicotiana_benthamiana",
        topologies: Sequence[str] = ("c_terminal", "n_terminal", "mir_loop"),
        n_variants_per_topology: int = 1,
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.48,
        seed: int = 42,
    ) -> list[dict[str, Any]]:
        """Generate 5DP VLP candidates and JSONL-ready yield-predictor records."""
        rng = random.Random(seed)
        results: list[dict[str, Any]] = []
        for topo in topologies:
            assembled = assemble_wssv_vlp_construct(hbcag_seq, vp28_seq, topology=topo)
            protein = assembled["protein_sequence"]
            for v_idx in range(max(1, n_variants_per_topology)):
                # Apply deterministic synonymous jitter on weights for diversity when N > 1
                if v_idx == 0:
                    w = dict(codon_weights)
                else:
                    w = {c: max(0.05, wt * rng.uniform(0.85, 1.15)) for c, wt in codon_weights.items()}
                opt = DPV240Optimizer(
                    forbidden_motifs=self.forbidden_motifs,
                    pause_positions=assembled["pause_positions"],
                    pause_target_adaptiveness=0.35,
                    pause_penalty_weight=3.5,
                ).optimize(
                    protein_sequence=protein,
                    codon_weights=w,
                    target_gc_min=target_gc_min,
                    target_gc_max=target_gc_max,
                    stop_codon="TAA",
                )
                cds = opt["sequence"]
                mfe = compute_5p_mfe_evidence(cds)
                results.append(
                    {
                        "candidate_id": f"wssv_{topo}_v{v_idx + 1}",
                        "host": host,
                        "topology": topo,
                        "protein_sequence": protein,
                        "cds_sequence": cds,
                        "predicted_cai": round(calculate_cai(cds, codon_weights), 4),
                        "gc_content_pct": round(calculate_gc(cds), 2),
                        "mrna_5prime_mfe_kcal_mol": mfe.get("mfe_kcal_mol") if isinstance(mfe, dict) else None,
                        "n_glycan_motif_count": len(assembled["n_glycan_motifs"]),
                        "n_glycan_motifs": assembled["n_glycan_motifs"],
                        "pause_positions": assembled["pause_positions"],
                        "has_pr1a": assembled["has_pr1a"],
                        "has_kdel": assembled["has_kdel"],
                        "wet_lab_yield_mg_kg": None,
                        "design_tier": "5DP_WSSV_VLP",
                        "epistemic_status": "UNVERIFIED_PUBLIC_CORPUS",
                        "biological_validation": "NOT_AVAILABLE",
                    }
                )
        return results

    @staticmethod
    def to_jsonl(candidates: Sequence[dict[str, Any]]) -> str:
        """Serialize candidate records to JSONL matching wssv_yield_prediction_candidates.jsonl."""
        return "\n".join(json.dumps(c, ensure_ascii=False) for c in candidates) + ("\n" if candidates else "")
