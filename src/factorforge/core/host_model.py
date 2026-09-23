"""Host Model Single Source of Truth for FactorForge Slate v2 (Job 293A).

Provides a unified abstraction for host codon frequencies, relative adaptiveness
weights (w_ij), reference GC bounds, and model provenance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.engines.profile.utils import get_data_path, load_codon_table

logger = logging.getLogger(__name__)

STOP_CODONS = {"TAA", "TAG", "TGA"}

HOST_METADATA_REGISTRY: Dict[str, Dict[str, Any]] = {
    "nbenthamiana": {
        "display_name": "Nicotiana benthamiana",
        "model_version": "NbeV1.1_HC_v3.7",
        "default_target_gc": 0.45,
        "reference_gc_range": (0.40, 0.47),
        "preferred_stop": "TAA",
    },
    "ntabacum": {
        "display_name": "Nicotiana tabacum (BY-2)",
        "model_version": "Nta_v1.0",
        "default_target_gc": 0.44,
        "reference_gc_range": (0.39, 0.46),
        "preferred_stop": "TAA",
    },
    "athaliana": {
        "display_name": "Arabidopsis thaliana",
        "model_version": "Ath_v1.0",
        "default_target_gc": 0.44,
        "reference_gc_range": (0.39, 0.46),
        "preferred_stop": "TAA",
    },
    "wolffia_globosa": {
        "display_name": "Wolffia globosa",
        "model_version": "Wgl_v1.0",
        "default_target_gc": 0.47,
        "reference_gc_range": (0.42, 0.50),
        "preferred_stop": "TAA",
    },
}


@dataclass
class HostModel:
    """Single source of truth for host codon usage, CAI weights, and provenance."""

    host_id: str
    display_name: str
    model_version: str
    codon_frequencies: Dict[str, float]
    codon_weights: Dict[str, float]
    aa_to_codons: Dict[str, List[str]]
    default_target_gc: float = 0.45
    reference_gc_range: Tuple[float, float] = (0.40, 0.47)
    extreme_gc_guard: Tuple[float, float] = (20.0, 80.0)
    preferred_stop: str = "TAA"
    provenance: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, host_id: str, data_dir: Optional[Path] = None) -> HostModel:
        """Load and construct a HostModel instance from data assets."""
        resolved_data_dir = data_dir or get_data_path()
        clean_host = str(host_id).strip().lower()

        try:
            raw_table = load_codon_table(clean_host, resolved_data_dir)
        except Exception as e:
            raise FileNotFoundError(f"Codon usage table for host '{clean_host}' not found: {e}")

        codons_dict = raw_table.get("codons", {})
        if not codons_dict or not isinstance(codons_dict, dict):
            raise ValueError(f"Invalid codon table format for host '{clean_host}'")

        # 1. Parse frequencies and build AA to codons map
        codon_frequencies: Dict[str, float] = {}
        aa_to_codons: Dict[str, List[str]] = {}

        for codon, entry in codons_dict.items():
            if not isinstance(codon, str) or len(codon) != 3:
                continue
            codon_upper = codon.upper()
            if isinstance(entry, dict):
                freq = float(entry.get("frequency", 0.0))
                aa = entry.get("aa") or STANDARD_GENETIC_CODE.get(codon_upper, "X")
            elif isinstance(entry, (int, float)):
                freq = float(entry)
                aa = STANDARD_GENETIC_CODE.get(codon_upper, "X")
            else:
                continue

            codon_frequencies[codon_upper] = freq
            if aa and aa != "*":
                aa_to_codons.setdefault(aa, []).append(codon_upper)

        # Ensure all standard amino acids are mapped
        for codon_upper, aa in STANDARD_GENETIC_CODE.items():
            if aa != "*" and codon_upper not in codon_frequencies:
                codon_frequencies[codon_upper] = 0.01
                if codon_upper not in aa_to_codons.get(aa, []):
                    aa_to_codons.setdefault(aa, []).append(codon_upper)

        # 2. Compute relative adaptiveness weights w_ij = f_ij / max_k(f_ik)
        codon_weights: Dict[str, float] = {}
        for aa, codons in aa_to_codons.items():
            max_freq = max([codon_frequencies.get(c, 0.0) for c in codons], default=0.01)
            for c in codons:
                f = codon_frequencies.get(c, 0.0)
                # Bound relative adaptiveness into [0.01, 1.0]
                codon_weights[c] = round(max(0.01, min(1.0, f / max_freq if max_freq > 0 else 0.5)), 4)

        # Metadata
        meta = HOST_METADATA_REGISTRY.get(clean_host, {
            "display_name": clean_host.capitalize(),
            "model_version": f"{clean_host}_v1.0",
            "default_target_gc": 0.45,
            "reference_gc_range": (0.35, 0.55),
            "preferred_stop": "TAA",
        })

        return cls(
            host_id=clean_host,
            display_name=meta["display_name"],
            model_version=meta["model_version"],
            codon_frequencies=codon_frequencies,
            codon_weights=codon_weights,
            aa_to_codons=aa_to_codons,
            default_target_gc=meta["default_target_gc"],
            reference_gc_range=meta["reference_gc_range"],
            extreme_gc_guard=(20.0, 80.0),
            preferred_stop=meta.get("preferred_stop", "TAA"),
            provenance={
                "source": raw_table.get("source", "FactorForge Data Registry"),
                "version": raw_table.get("version", meta["model_version"]),
            },
        )

    @classmethod
    def is_host_available(cls, host_id: str, data_dir: Optional[Path] = None) -> bool:
        """Check if a host codon model is available in the repository without throwing."""
        try:
            cls.load(host_id, data_dir)
            return True
        except Exception:
            return False
