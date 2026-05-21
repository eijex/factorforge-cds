"""
Multidimensional Scoring for FactorForge v2.
Composite scoring function with optional ViennaRNA MFE integration.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

# Optimal GC range for N. benthamiana
GC_OPT_MIN = 41.0
GC_OPT_MAX = 44.0
GC_OPT_MID = 42.5

# ViennaRNA availability cache
_vienna_available: bool | None = None


@dataclass
class ScoringConfig:
    """Scoring weight configuration for composite optimization score."""

    w_cai: float = 0.5
    w_gc: float = 0.3
    w_mfe: float = 0.2
    w_dinuc: float = 0.0  # CpG/TpA dinucleotide penalty (opt-in, default off)
    gc_opt: float = GC_OPT_MID
    use_mfe: bool = True

    def __post_init__(self) -> None:
        """Normalize weights to sum to 1.0."""
        self._normalize()

    def _normalize(self) -> None:
        """Ensure active weights sum to 1.0."""
        total = (
            self.w_cai
            + self.w_gc
            + (self.w_mfe if self.use_mfe else 0.0)
            + self.w_dinuc
        )
        if total > 0:
            self.w_cai /= total
            self.w_gc /= total
            if self.use_mfe:
                self.w_mfe /= total
            else:
                self.w_mfe = 0.0
            self.w_dinuc /= total


# Pre-defined scoring configs per optimization profile
PROFILE_SCORING_CONFIGS: dict[str, ScoringConfig] = {
    "balanced": ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, gc_opt=GC_OPT_MID),
    "high_cai": ScoringConfig(w_cai=0.8, w_gc=0.1, w_mfe=0.1, gc_opt=GC_OPT_MID),
    "gc_target": ScoringConfig(w_cai=0.1, w_gc=0.7, w_mfe=0.2, gc_opt=50.0),
    "assembly_friendly": ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2, gc_opt=GC_OPT_MID),
    "ramp": ScoringConfig(w_cai=0.4, w_gc=0.3, w_mfe=0.3, gc_opt=GC_OPT_MID),
    # TRV viral-delivery profile — Li et al. (2026): prioritize MFE and viral-context GC target.
    "viral_delivery": ScoringConfig(w_cai=0.35, w_gc=0.25, w_mfe=0.40, gc_opt=47.5, use_mfe=True),
}


def _check_vienna_available() -> bool:
    """Check if ViennaRNA Python bindings are available."""
    global _vienna_available
    if _vienna_available is None:
        try:
            import RNA  # noqa: F401

            _vienna_available = True
            logger.debug("ViennaRNA Python bindings available")
        except ImportError:
            _vienna_available = False
            logger.debug("ViennaRNA not available; MFE scoring disabled")
    return _vienna_available


def calculate_mfe(sequence: str) -> float | None:
    """
    Calculate minimum free energy (MFE) using ViennaRNA.

    Args:
        sequence: DNA or RNA sequence.

    Returns:
        MFE in kcal/mol, or None if ViennaRNA is not available.
    """
    if not _check_vienna_available():
        return None

    try:
        import RNA

        # Convert DNA to RNA (T → U)
        rna_seq = sequence.upper().replace("T", "U")
        _, mfe = RNA.fold(rna_seq)
        return float(mfe)
    except Exception as exc:
        logger.debug(f"MFE calculation failed: {exc}")
        return None


def normalize_mfe(mfe: float, seq_length: int) -> float:
    """
    Normalize MFE to 0-1 range where 1 = no structure (favorable).

    Uses empirical scaling: MFE per nucleotide typically ranges from
    -0.5 to 0.0 kcal/mol/nt for mRNA coding sequences.

    Args:
        mfe: Minimum free energy in kcal/mol.
        seq_length: Sequence length in nucleotides.

    Returns:
        Normalized MFE score (0-1, higher = less structured = better for translation).
    """
    if seq_length == 0:
        return 0.5

    mfe_per_nt = mfe / seq_length
    # Clamp to expected range [-0.5, 0.0]
    clamped = max(-0.5, min(0.0, mfe_per_nt))
    # Map to [0, 1] where 0.0 kcal/mol/nt → 1.0 and -0.5 → 0.0
    return 1.0 + (clamped / 0.5)


def calculate_dinucleotide_score(sequence: str) -> float:
    """Calculate a dinucleotide avoidance score (0-1, higher = fewer CpG/TpA).

    Combines CpG and TpA observed/expected ratios. A sequence with no CpG
    and no TpA scores 1.0; high density scores toward 0.0.

    Args:
        sequence: DNA sequence.

    Returns:
        Dinucleotide avoidance score (0-1).
    """
    from factorforge.engines.v2.utils import calculate_dinucleotide_ratio

    if len(sequence) < 6:
        return 1.0

    cpg_ratio = calculate_dinucleotide_ratio(sequence, "CG")
    tpa_ratio = calculate_dinucleotide_ratio(sequence, "TA")

    # Score: 1.0 when ratio=0, 0.0 when ratio>=2.0
    cpg_score = max(0.0, 1.0 - cpg_ratio / 2.0)
    tpa_score = max(0.0, 1.0 - tpa_ratio / 2.0)

    return (cpg_score + tpa_score) / 2.0


def calculate_composite_score(
    cai: float,
    gc: float,
    sequence: str | None = None,
    config: ScoringConfig | None = None,
    profile: str | None = None,
    **kwargs: Any,
) -> float:
    """
    Calculate multidimensional composite score.

    S = w1*CAI + w2*(1 - |GC - GC_opt|/50) + w3*MFE_norm + w4*dinuc_score

    Args:
        cai: Codon Adaptation Index (0-1).
        gc: GC content percentage (0-100).
        sequence: DNA sequence for optional MFE calculation.
        config: Explicit ScoringConfig. Overrides profile.
        profile: Profile name for preset config lookup.
        **kwargs: Additional parameters (e.g., target_gc for gc_target profile).

    Returns:
        Composite score (0-1).
    """
    # Resolve config
    if config is None:
        profile_name = (profile or "balanced").lower()
        config = PROFILE_SCORING_CONFIGS.get(profile_name)
        if config is None:
            config = PROFILE_SCORING_CONFIGS["balanced"]

    # Allow target_gc override for gc_target profile
    gc_opt = float(kwargs.get("target_gc", config.gc_opt))

    # Component 1: CAI (already 0-1)
    cai_score = max(0.0, min(1.0, cai))

    # Component 2: GC proximity to optimum
    gc_score = max(0.0, 1.0 - abs(gc - gc_opt) / 50.0)

    # Component 3: MFE (optional)
    mfe_score = 0.5  # neutral default
    actual_w_mfe = config.w_mfe

    if config.use_mfe and sequence is not None and _check_vienna_available():
        mfe = calculate_mfe(sequence)
        if mfe is not None:
            mfe_score = normalize_mfe(mfe, len(sequence))
        else:
            actual_w_mfe = 0.0
    else:
        actual_w_mfe = 0.0

    # Component 4: Dinucleotide avoidance (opt-in, default weight 0.0)
    dinuc_score = 0.5  # neutral default
    actual_w_dinuc = config.w_dinuc
    if actual_w_dinuc > 0 and sequence is not None:
        dinuc_score = calculate_dinucleotide_score(sequence)
    elif actual_w_dinuc > 0:
        actual_w_dinuc = 0.0  # Cannot compute without sequence

    # Compute weighted score (re-normalize if MFE/dinuc disabled)
    w_total = config.w_cai + config.w_gc + actual_w_mfe + actual_w_dinuc
    if w_total == 0:
        return 0.0

    score = (
        (config.w_cai / w_total) * cai_score
        + (config.w_gc / w_total) * gc_score
        + (actual_w_mfe / w_total) * mfe_score
        + (actual_w_dinuc / w_total) * dinuc_score
    )

    return round(score, 3)
