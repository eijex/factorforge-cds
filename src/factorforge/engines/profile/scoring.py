"""
Multidimensional Scoring for FactorForge profile engine.
Composite scoring function with optional ViennaRNA MFE integration.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

# GC band for N. benthamiana codon-optimized sequences.
# Benchmark (internal, n=49): balanced profile output average GC% = 60.1%
# (range 55-71%). The genome-wide average (~42%) reflects all genes, not the
# high-expression codon table which exhibits 3rd-position GC bias.
# These constants define the acceptable band — sequences within [GC_OPT_MIN, GC_OPT_MAX]
# receive full GC score; outside the band the score decays linearly.
GC_OPT_MIN = 55.0
GC_OPT_MAX = 65.0
GC_OPT_MID = 60.0  # kept for gc_target point-scoring and viral_delivery centering
GC_DECAY_WIDTH = 20.0  # percentage points outside band before score reaches 0.0

# ViennaRNA availability cache
_vienna_available: bool | None = None


@dataclass
class ScoringConfig:
    """Scoring weight configuration for composite optimization score."""

    w_cai: float = 0.5
    w_gc: float = 0.3
    w_mfe: float = 0.2
    w_dinuc: float = 0.0  # CpG/TpA dinucleotide penalty (opt-in, default off)
    cpg_weight: float = 0.0  # plant default: CpG inactive
    tpa_weight: float = 1.0  # plant default: TpA active
    w_syncodonlm: float = 0.0  # SynCodonLM quality score (opt-in, default off)
    gc_opt: float = GC_OPT_MID  # no longer used by calculate_composite_score (superseded by
                                # gc_min/gc_max band); retained for external API compatibility
    gc_min: float = GC_OPT_MIN  # acceptable band lower boundary
    gc_max: float = GC_OPT_MAX  # acceptable band upper boundary
    gc_decay_width: float = GC_DECAY_WIDTH  # % points outside band before score → 0
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
            + self.w_syncodonlm
        )
        if total > 0:
            self.w_cai /= total
            self.w_gc /= total
            if self.use_mfe:
                self.w_mfe /= total
            else:
                self.w_mfe = 0.0
            self.w_dinuc /= total
            self.w_syncodonlm /= total


# Pre-defined scoring configs per optimization profile
PROFILE_SCORING_CONFIGS: dict[str, ScoringConfig] = {
    "balanced": ScoringConfig(w_cai=0.5, w_gc=0.3, w_mfe=0.2),
    # high_cai: CAI 1.0 mimics naturally high-expression N. benthamiana genes (golden_set).
    # For very long or structurally complex proteins under extreme agroinfiltration overexpression,
    # consider the ramp profile instead (N-terminal codon deoptimization reduces ribosome stalling).
    "high_cai": ScoringConfig(w_cai=0.8, w_gc=0.1, w_mfe=0.1),
    # gc_target: gc_min/gc_max are overridden dynamically from target_gc kwarg in
    # calculate_composite_score — the config defaults here are not used for band scoring.
    "gc_target": ScoringConfig(w_cai=0.1, w_gc=0.7, w_mfe=0.2),
    # assembly_friendly: CAI pressure reduced vs balanced; GC/MFE weights raised to
    # reflect Type IIS restriction-site avoidance priority (Golden Gate / MoClo).
    # w_gc scores GC band compliance (global GC%), NOT local GC uniformity.
    # Window-level GC variance and repeat-pattern penalties are not yet implemented.
    "assembly_friendly": ScoringConfig(w_cai=0.3, w_gc=0.4, w_mfe=0.3),
    "ramp": ScoringConfig(w_cai=0.4, w_gc=0.3, w_mfe=0.3),
    # TRV viral-delivery profile — GC band centered on TRV genome composition (~47.5%).
    # MFE weighted at 0.30 (Peccoud et al. 2024, PMC11718241).
    "viral_delivery": ScoringConfig(
        w_cai=0.35, w_gc=0.35, w_mfe=0.30,
        gc_opt=47.5, gc_min=37.5, gc_max=57.5,
        use_mfe=True,
    ),
    "ml_enhanced": ScoringConfig(w_cai=0.35, w_gc=0.25, w_mfe=0.15, w_syncodonlm=0.25),
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


def gc_band_score(
    gc: float,
    gc_min: float,
    gc_max: float,
    decay_width: float = GC_DECAY_WIDTH,
) -> float:
    """Score GC content against an acceptable band.

    Returns 1.0 inside [gc_min, gc_max]; linearly decays to 0.0 after
    decay_width percentage points outside the band.

    Args:
        gc: GC content percentage (0-100).
        gc_min: Lower boundary of acceptable band.
        gc_max: Upper boundary of acceptable band.
        decay_width: Percentage points outside band before score reaches 0.0.

    Examples:
        gc_min=55, gc_max=65, decay_width=20:
          gc=60 → 1.00  (inside band)
          gc=70 → 0.75  (5 pts above gc_max)
          gc=80 → 0.25  (15 pts above gc_max)
          gc=85 → 0.00  (20 pts above gc_max)
    """
    if gc_min <= gc <= gc_max:
        return 1.0
    distance = (gc_min - gc) if gc < gc_min else (gc - gc_max)
    return max(0.0, 1.0 - distance / decay_width)


def calculate_dinucleotide_score(
    sequence: str,
    cpg_weight: float = 0.0,
    tpa_weight: float = 1.0,
) -> float:
    """Score dinucleotide avoidance.

    Plant default: CpG inactive (cpg_weight=0.0), only TpA is penalized.
    Mammalian opt-in: set cpg_weight=1.0 and tpa_weight=1.0 to penalize both.
    """
    from factorforge.analysis.metrics import calculate_dinucleotide_score as _score

    return _score(sequence, cpg_weight=cpg_weight, tpa_weight=tpa_weight)


def calculate_composite_score(
    cai: float,
    gc: float,
    sequence: str | None = None,
    config: ScoringConfig | None = None,
    profile: str | None = None,
    **kwargs: Any,
) -> float:
    """Calculate multidimensional composite score.

    S = w1*CAI + w2*gc_band_score + w3*MFE_norm
        + w4*dinuc_score + w5*SynCodonLM_score

    GC scoring uses a band function: sequences inside [gc_min, gc_max] receive
    full score (1.0); outside the band the score decays linearly to 0.0 over
    gc_decay_width percentage points. For gc_target profile, the band is
    centred on the caller-supplied target_gc (±5 pp).

    Args:
        cai: Codon Adaptation Index (0-1).
        gc: GC content percentage (0-100).
        sequence: DNA sequence for optional MFE, dinucleotide, and SynCodonLM calculation.
        config: Explicit ScoringConfig. Overrides profile.
        profile: Profile name for preset config lookup.
        **kwargs: target_gc (float) — point target for gc_target profile.

    Returns:
        Composite score (0-1).
    """
    # Resolve config
    if config is None:
        profile_name = (profile or "balanced").lower()
        config = PROFILE_SCORING_CONFIGS.get(profile_name)
        if config is None:
            config = PROFILE_SCORING_CONFIGS["balanced"]

    # Component 1: CAI (already 0-1)
    cai_score = max(0.0, min(1.0, cai))

    # Component 2: GC band scoring
    # gc_target profile: caller supplies target_gc; use a ±5 pp band around it.
    # All other profiles: use the band defined in ScoringConfig (gc_min/gc_max).
    if "target_gc" in kwargs:
        tgt = float(kwargs["target_gc"])
        gc_score = gc_band_score(gc, tgt - 5.0, tgt + 5.0, config.gc_decay_width)
    else:
        gc_score = gc_band_score(gc, config.gc_min, config.gc_max, config.gc_decay_width)

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
        dinuc_score = calculate_dinucleotide_score(
            sequence,
            cpg_weight=config.cpg_weight,
            tpa_weight=config.tpa_weight,
        )
    elif actual_w_dinuc > 0:
        actual_w_dinuc = 0.0  # Cannot compute without sequence

    # Component 5: SynCodonLM score (opt-in, default weight 0.0)
    syncodonlm_score = 0.5  # neutral default
    actual_w_syncodonlm = config.w_syncodonlm
    if actual_w_syncodonlm > 0 and sequence is not None:
        from factorforge.engines.profile.scoring_ml import calculate_syncodonlm_score

        organism = str(kwargs.get("organism", "Nicotiana_benthamiana"))
        syncodonlm_score = calculate_syncodonlm_score(sequence, organism=organism)
    elif actual_w_syncodonlm > 0:
        actual_w_syncodonlm = 0.0  # Cannot compute without sequence

    # Compute weighted score (re-normalize if optional components are disabled)
    w_total = config.w_cai + config.w_gc + actual_w_mfe + actual_w_dinuc + actual_w_syncodonlm
    if w_total == 0:
        return 0.0

    score = (
        (config.w_cai / w_total) * cai_score
        + (config.w_gc / w_total) * gc_score
        + (actual_w_mfe / w_total) * mfe_score
        + (actual_w_dinuc / w_total) * dinuc_score
        + (actual_w_syncodonlm / w_total) * syncodonlm_score
    )

    return round(score, 3)


def compute_mfe_evidence(
    sequence: str | None,
    config: ScoringConfig | None = None,
    profile: str | None = None,
) -> dict[str, Any]:
    """Return MFE provenance metadata for a scored sequence.

    Mirrors the MFE branch of ``calculate_composite_score`` WITHOUT changing the
    score. Its purpose is honesty: when MFE is not computed (e.g. ViennaRNA is
    unavailable in the deployment, as on Vercel), callers must be able to tell
    that ``mfe_kcal_mol`` is absent rather than a genuine 0.0.

    Returns a dict with:
        mfe_kcal_mol: float | None   (None when not computed)
        mfe_status:   "computed" | "not_computed"
        mfe_used:     bool           (whether MFE contributed to the score)
        mfe_warning:  str | None     (reason when not used)
        score_components: {cai_used, gc_used, mfe_used}
    """
    if config is None:
        profile_name = (profile or "balanced").lower()
        config = PROFILE_SCORING_CONFIGS.get(profile_name) or PROFILE_SCORING_CONFIGS["balanced"]

    mfe_value: float | None = None
    reason: str | None = None

    if not config.use_mfe:
        reason = "MFE scoring is disabled for this profile."
    elif sequence is None:
        reason = "MFE was not computed because no sequence was provided."
    elif not _check_vienna_available():
        reason = "MFE was not computed because ViennaRNA is unavailable in this environment."
    else:
        mfe_value = calculate_mfe(sequence)
        if mfe_value is None:
            reason = "MFE computation failed for this sequence."

    mfe_used = mfe_value is not None
    return {
        "mfe_kcal_mol": round(mfe_value, 2) if mfe_used else None,
        "mfe_status": "computed" if mfe_used else "not_computed",
        "mfe_used": mfe_used,
        "mfe_warning": None if mfe_used else reason,
        "score_components": {"cai_used": True, "gc_used": True, "mfe_used": mfe_used},
    }
