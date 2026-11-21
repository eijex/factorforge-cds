"""Experimental 5DP Prototype: Multi-Part Fusion Assembly & Plant PTM Layer.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Implements a research-preview 5DP scorer/generator for multi-domain fusion
and chimeric protein constructs, tracking assembly junction pause windows and
canonical N-glycosylation sequons (N-X-S/T, X != P).
"""

from __future__ import annotations

import json
import random
from collections.abc import Sequence
from typing import Any

from factorforge.analysis.metrics import calculate_cai, calculate_gc
from factorforge.engines.dp_v2_4 import DPV240Optimizer
from factorforge.engines.profile.scoring import compute_5p_mfe_evidence


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


def assemble_fusion_construct(
    scaffold_seq: str,
    cargo_seq: str,
    topology: str = "c_terminal",
    signal_peptide: str = "",
    linker: str = "GGGGS",
    c_tag: str = "",
) -> dict[str, Any]:
    """Assemble a multi-part fusion protein and compute junction pause indices."""
    scaf = "".join(scaffold_seq.upper().split()).rstrip("*")
    cargo = "".join(cargo_seq.upper().split()).rstrip("*")
    if not scaf or not cargo:
        raise ValueError("Both scaffold_seq and cargo_seq must be non-empty.")

    prefix = "".join(signal_peptide.upper().split()).rstrip("*")
    suffix = "".join(c_tag.upper().split()).rstrip("*")

    topo = topology.strip().lower()
    if topo == "c_terminal":
        core = scaf + linker + cargo
        boundary_aa = len(prefix) + len(scaf)
        steric_note = "C-terminal fusion; flexible linker separates scaffold from cargo."
    elif topo == "n_terminal":
        core = cargo + linker + scaf
        boundary_aa = len(prefix) + len(cargo)
        steric_note = "N-terminal fusion; cargo precedes scaffold assembly domain."
    elif topo == "internal_loop":
        mid = len(scaf) // 2
        core = scaf[:mid] + linker + cargo + linker + scaf[mid:]
        boundary_aa = len(prefix) + mid
        steric_note = "Internal-loop insertion; surface-exposed display (steric check required)."
    else:
        raise ValueError("topology must be 'c_terminal', 'n_terminal', or 'internal_loop'.")

    full_protein = prefix + core + suffix
    # 5DP junction boundary + exit-tunnel offset (+35 aa) pause window
    pause_center = min(len(full_protein), boundary_aa + 35)
    pause_positions = [p for p in (pause_center - 1, pause_center, pause_center + 1) if 1 <= p <= len(full_protein)]

    return {
        "topology": topo,
        "protein_sequence": full_protein,
        "length_aa": len(full_protein),
        "fusion_boundary_aa": boundary_aa,
        "pause_positions": pause_positions,
        "n_glycan_motifs": detect_n_glycosylation_motifs(cargo),
        "has_signal_peptide": bool(prefix),
        "has_c_tag": bool(suffix),
        "steric_note": steric_note,
    }


class Assembly5DPOptimizer:
    """Experimental 5DP candidate generator for multi-part fusion constructs."""

    def __init__(self, forbidden_motifs: list[str] | None = None) -> None:
        self.forbidden_motifs = forbidden_motifs or ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"]

    def generate_candidates(
        self,
        scaffold_seq: str,
        cargo_seq: str,
        codon_weights: dict[str, float],
        host: str = "nicotiana_benthamiana",
        topologies: Sequence[str] = ("c_terminal", "n_terminal", "internal_loop"),
        signal_peptide: str = "",
        linker: str = "GGGGS",
        c_tag: str = "",
        n_variants_per_topology: int = 1,
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.48,
        seed: int = 42,
    ) -> list[dict[str, Any]]:
        """Generate 5DP assembly candidates and JSONL-ready candidate records."""
        rng = random.Random(seed)
        results: list[dict[str, Any]] = []
        for topo in topologies:
            assembled = assemble_fusion_construct(
                scaffold_seq,
                cargo_seq,
                topology=topo,
                signal_peptide=signal_peptide,
                linker=linker,
                c_tag=c_tag,
            )
            protein = assembled["protein_sequence"]
            for v_idx in range(max(1, n_variants_per_topology)):
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
                        "candidate_id": f"assembly_{topo}_v{v_idx + 1}",
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
                        "has_signal_peptide": assembled["has_signal_peptide"],
                        "has_c_tag": assembled["has_c_tag"],
                        "wet_lab_yield_mg_kg": None,
                        "design_tier": "5DP_ASSEMBLY",
                        "epistemic_status": "UNVERIFIED_PUBLIC_CORPUS",
                        "biological_validation": "NOT_AVAILABLE",
                    }
                )
        return results

    @staticmethod
    def to_jsonl(candidates: Sequence[dict[str, Any]]) -> str:
        """Serialize candidate records to JSONL."""
        return "\n".join(json.dumps(c, ensure_ascii=False) for c in candidates) + ("\n" if candidates else "")
