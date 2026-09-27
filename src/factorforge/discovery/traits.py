"""Multi-Dimensional Trait Vector Extraction for FactorForge Discovery Slate (Job 283A).

Extracts biophysical, adaptation, structural, and synthesis traits for candidate sequences.
Evidence Tier: COMPUTATIONAL_HEURISTIC.
"""

from __future__ import annotations
from typing import Optional

from factorforge.analysis.metrics import (
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    detect_repeats,
    load_codon_usage_table,
)
from factorforge.discovery.schemas import TraitVector
from factorforge.engines.profile.scoring import (
    calculate_mfe,
)


class TraitVectorExtractor:
    """Computes multi-dimensional trait vectors for candidate CDS sequences."""

    def __init__(self, host: str = "nbenthamiana") -> None:
        self.host = host
        self.codon_table = load_codon_usage_table()

    def extract(self, sequence_dna: str) -> TraitVector:
        """Extracts the 8-dimensional trait vector from a sequence."""
        seq = sequence_dna.strip().upper()

        # 1. Golden Set CAI
        cai_val = calculate_cai(seq, self.codon_table.codon_weights)

        # 2. Global GC%
        gc_val = calculate_gc(seq)

        # 3. 5' Initiation MFE (kcal/mol)
        mfe_val: Optional[float] = None
        if len(seq) >= 30:
            initiation_window = seq[: min(len(seq), 60)]
            mfe_val = calculate_mfe(initiation_window)

        # 4. Local 50-bp GC min and max
        if len(seq) >= 50:
            windows = calculate_gc_windows(seq, window_size=50)
            gc_vals = [float(w["gc"]) for w in windows]
            local_gc_min = min(gc_vals) if gc_vals else gc_val
            local_gc_max = max(gc_vals) if gc_vals else gc_val
        else:
            local_gc_min = gc_val
            local_gc_max = gc_val

        # 5. Maximum Homopolymer Run Length
        max_homopolymer = 1
        current_run = 1
        for i in range(1, len(seq)):
            if seq[i] == seq[i - 1]:
                current_run += 1
                if current_run > max_homopolymer:
                    max_homopolymer = current_run
            else:
                current_run = 1

        # 6. Synthesis Penalty (0.0 to 1.0, lower is better)
        repeats = detect_repeats(seq)
        repeat_penalty = min(0.5, len(repeats) * 0.1)

        # Local GC extreme penalty
        gc_extreme_penalty = 0.0
        if local_gc_min < 20.0:
            gc_extreme_penalty += (20.0 - local_gc_min) * 0.01
        if local_gc_max > 70.0:
            gc_extreme_penalty += (local_gc_max - 70.0) * 0.01

        synthesis_penalty = min(1.0, round(repeat_penalty + gc_extreme_penalty, 4))

        # 7. Codon-Context Prior Score (0.0 to 1.0 heuristic)
        # Evaluates harmony between sequence codon distribution and reference weights
        prior_score = round(
            max(0.0, min(1.0, (cai_val * 0.7) + (1.0 - synthesis_penalty) * 0.3)), 4
        )

        return TraitVector(
            cai_golden_set=round(cai_val, 4),
            global_gc_percent=round(gc_val, 2),
            initiation_mfe_kcal_mol=round(mfe_val, 2) if mfe_val is not None else None,
            local_50bp_gc_min=round(local_gc_min, 1),
            local_50bp_gc_max=round(local_gc_max, 1),
            homopolymer_max_run=max_homopolymer,
            synthesis_penalty=synthesis_penalty,
            codon_context_prior=prior_score,
        )
