"""Optional, validated initiation-window observation with per-call model settings."""

from __future__ import annotations
from dataclasses import asdict, dataclass, field
import math
import time
from typing import Any, Dict, Optional

DEFAULT_INITIATION_WINDOW_NT = 48
DEFAULT_TEMPERATURE_C = 37.0
# Retain the existing software advisory; it is not a calibrated plant cutoff.
ADVISORY_MFE_THRESHOLD_KCAL_MOL = -3.0
TRUTHFULNESS_NOTE = (
    "Thermodynamic MFE is an in-silico RNA secondary structure proxy for the 5' initiation region. "
    "It does not guarantee translation velocity, ribosomal initiation frequency, or in-vivo yield."
)


@dataclass
class InitiationMFEResult:
    status: str
    evidence_class: str = "MODEL_PREDICTION"
    engine: str = "ViennaRNA_RNAfold"
    input_scope: str = "cds_5p_only"
    window_start: int = 0
    window_end: int = 0
    folded_sequence_length: int = 0
    temperature_c: float = DEFAULT_TEMPERATURE_C
    mfe_kcal_mol: Optional[float] = None
    secondary_structure: Optional[str] = None
    runtime_ms: float = 0.0
    reason: Optional[str] = None
    advisory_flag: Optional[str] = None
    truthfulness_note: str = TRUTHFULNESS_NOTE
    parameter_set: str = "engine_current_parameters"
    model_details: Dict[str, Any] = field(default_factory=dict)
    advisory_threshold: Optional[float] = ADVISORY_MFE_THRESHOLD_KCAL_MOL
    upstream_bases_used: int = 0
    coordinate_system: str = "zero_based_relative_to_cds_start_end_exclusive"

    def to_dict(self):
        return asdict(self)


def evaluate_initiation_mfe(
    cds: str,
    upstream_context: Optional[str] = None,
    window_nt: int = DEFAULT_INITIATION_WINDOW_NT,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    advisory_threshold: Optional[float] = ADVISORY_MFE_THRESHOLD_KCAL_MOL,
) -> InitiationMFEResult:
    started = time.perf_counter()
    result = InitiationMFEResult(status="ERROR")

    def finish(status, reason=None):
        result.status = status
        result.reason = reason
        result.runtime_ms = round((time.perf_counter() - started) * 1000, 2)
        return result

    if type(window_nt) is not int or not 1 <= window_nt <= 512:
        return finish("ERROR", "window_nt must be an integer in [1, 512]")
    if (
        isinstance(temperature_c, bool)
        or not isinstance(temperature_c, (int, float))
        or not math.isfinite(temperature_c)
        or not 0 <= temperature_c <= 100
    ):
        return finish("ERROR", "temperature must be finite and in [0, 100] Celsius")
    if advisory_threshold is not None and (
        isinstance(advisory_threshold, bool)
        or not isinstance(advisory_threshold, (int, float))
        or not math.isfinite(advisory_threshold)
    ):
        return finish("ERROR", "advisory threshold must be finite or null")
    if not isinstance(cds, str) or (
        upstream_context is not None and not isinstance(upstream_context, str)
    ):
        return finish("ERROR", "CDS and upstream context must be strings")
    clean = "".join(cds.upper().split())
    up = "".join((upstream_context or "").upper().split())
    if not clean or set(clean + up) - set("ACGTU"):
        return finish("ERROR", "Folding requires nonempty, unambiguous nucleotides")
    # Use up to four supplied bases directly adjacent to CDS start, never far-upstream-only.
    used = up[-min(4, window_nt - 1) :] if up and window_nt > 1 else ""
    fold_window = (used + clean[: window_nt - len(used)]).replace("T", "U")
    result.input_scope = "upstream_flank_plus_cds_5p" if used else "cds_5p_only"
    result.window_start = -len(used)
    result.window_end = len(fold_window) - len(used)
    result.upstream_bases_used = len(used)
    result.folded_sequence_length = len(fold_window)
    result.temperature_c = float(temperature_c)
    result.advisory_threshold = advisory_threshold
    try:
        import RNA
    except ImportError:
        result.advisory_flag = "UNMEASURED"
        return finish(
            "NOT_AVAILABLE", "ViennaRNA package is not installed in the current environment"
        )
    try:
        result.engine = f"ViennaRNA_{getattr(RNA, '__version__', 'unknown')}"
        if hasattr(RNA, "md") and hasattr(RNA, "fold_compound"):
            md = RNA.md()
            md.temperature = float(temperature_c)
            result.model_details = {
                k: getattr(md, k)
                for k in ("temperature", "dangles", "noLP", "maxBPspan")
                if hasattr(md, k)
            }
            structure, energy = RNA.fold_compound(fold_window, md).mfe()
        else:
            result.advisory_flag = "UNMEASURED"
            return finish("NOT_AVAILABLE", "RNA engine lacks controlled per-call model settings")
        if (
            isinstance(energy, bool)
            or not math.isfinite(float(energy))
            or not isinstance(structure, str)
            or len(structure) != len(fold_window)
            or set(structure) - set(".()")
        ):
            raise ValueError("invalid engine result")
        balance = 0
        for char in structure:
            balance += (char == "(") - (char == ")")
            if balance < 0:
                raise ValueError("invalid structure")
        if balance:
            raise ValueError("invalid structure")
        result.mfe_kcal_mol = round(float(energy), 2)
        result.secondary_structure = structure
        if advisory_threshold is not None and float(energy) < advisory_threshold:
            result.advisory_flag = (
                "Below legacy software advisory threshold; plant efficacy cutoff is uncalibrated"
            )
        return finish("AVAILABLE")
    except Exception:
        result.mfe_kcal_mol = None
        result.secondary_structure = None
        result.advisory_flag = "ERROR"
        return finish("ERROR", "RNA folding failed or returned invalid observations")
