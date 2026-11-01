"""ESM-2 Neural Contact & Flexibility Propensity Recommendation Adapter.

Provides fast, sequence-only structural propensity profiling using pre-trained
ESM-2 representations (esm2_t6_8M_UR50D).

Epistemic & Governance Constraints:
- Biological validation: NOT_AVAILABLE.
- Classification: UNVERIFIED_PUBLIC_CORPUS.
- Computational contact density valleys do not prove ribosomal pausing benefit or folding kinetics.
- Strict human-in-the-loop invariant: reports always initialize with reviewed=False.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from typing import Any, Dict, List, Optional

from factorforge.structure.auto_annotator import CandidateRegion

ESM_ANNOTATOR_VERSION = "1.0.0"
ESM_LIMITATIONS_STATEMENT = (
    "Computational structural contact density and representation profiles from ESM-2 (esm2_t6_8M_UR50D) "
    "are exploratory predictions based on public pre-trained neural networks. "
    "Valleys in contact density or token entropy do not establish confirmed biological loops, "
    "ribosomal dwell kinetics, or in vivo folding efficiency. "
    "Biological validation: NOT_AVAILABLE. Classification: UNVERIFIED_PUBLIC_CORPUS. "
    "Explicit human review is required before use in CDS design."
)


@dataclass
class ESMRecommendationReport:
    """Full recommendation payload from ESM-2 contact density analysis."""

    protein_length: int
    model_name: str
    contact_density_min: float
    contact_density_mean: float
    contact_density_max: float
    candidate_regions: List[CandidateRegion] = field(default_factory=list)
    recommended_positions: List[int] = field(default_factory=list)
    coverage: float = 0.0
    tool_version: str = ESM_ANNOTATOR_VERSION
    limitations: str = ESM_LIMITATIONS_STATEMENT
    status: str = "AVAILABLE"
    reason: Optional[str] = None
    provenance_basis: str = "UNVERIFIED_PUBLIC_CORPUS"
    biological_validation: str = "NOT_AVAILABLE"
    reviewed: bool = False  # Hard invariant: Always False until explicit human confirmation

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class ESMStructureAnnotator:
    """Lightweight sequence-only structural propensity annotator powered by ESM-2."""

    def __init__(self, model_name: str = "esm2_t6_8M_UR50D") -> None:
        self.model_name = model_name
        self._model = None
        self._alphabet = None
        self._batch_converter = None
        self._load_attempted = False

    def _ensure_model_loaded(self) -> bool:
        if self._model is not None:
            return True
        if self._load_attempted:
            return False

        self._load_attempted = True
        try:
            import esm  # type: ignore

            loader = getattr(esm.pretrained, self.model_name, None)
            if loader is None:
                return False
            model, alphabet = loader()
            model.eval()
            self._model = model
            self._alphabet = alphabet
            self._batch_converter = alphabet.get_batch_converter()
            return True
        except Exception:
            return False

    def is_available(self) -> bool:
        """Check if ESM model can be imported and loaded."""
        return self._ensure_model_loaded()

    def annotate(
        self,
        protein_sequence: str,
        threshold_std: float = 0.75,
        min_region_len: int = 3,
        smoothing_window: int = 5,
    ) -> ESMRecommendationReport:
        """Analyze protein sequence with ESM-2 to identify contact density valleys."""
        clean_prot = "".join(protein_sequence.upper().split()).rstrip("*")
        n_res = len(clean_prot)

        if n_res == 0:
            return ESMRecommendationReport(
                protein_length=0,
                model_name=self.model_name,
                contact_density_min=0.0,
                contact_density_mean=0.0,
                contact_density_max=0.0,
                status="ERROR",
                reason="Protein sequence is empty",
            )

        if not self._ensure_model_loaded():
            return ESMRecommendationReport(
                protein_length=n_res,
                model_name=self.model_name,
                contact_density_min=0.0,
                contact_density_mean=0.0,
                contact_density_max=0.0,
                status="UNAVAILABLE",
                reason=f"ESM model '{self.model_name}' could not be loaded or fair-esm is not installed.",
            )

        try:
            import torch  # type: ignore

            data = [("query", clean_prot)]
            _, _, batch_tokens = self._batch_converter(data)

            with torch.no_grad():
                results = self._model(batch_tokens, repr_layers=[6], return_contacts=True)

            contacts = results["contacts"][0].cpu().numpy()  # shape (L, L)

            # Compute non-local contact density (|i - j| >= 6)
            non_local = [0.0] * n_res
            for i in range(n_res):
                s = 0.0
                for j in range(n_res):
                    if abs(i - j) >= 6:
                        s += float(contacts[i, j])
                non_local[i] = s

            # Moving average smoothing
            half_w = smoothing_window // 2
            smoothed = [0.0] * n_res
            for i in range(n_res):
                w_start = max(0, i - half_w)
                w_end = min(n_res, i + half_w + 1)
                sub = non_local[w_start:w_end]
                smoothed[i] = sum(sub) / len(sub)

            mean_c = sum(smoothed) / n_res
            var_c = sum((x - mean_c) ** 2 for x in smoothed) / max(1, n_res)
            std_c = math.sqrt(var_c)

            c_min = min(smoothed)
            c_max = max(smoothed)
            valley_threshold = mean_c - (threshold_std * std_c)

            # Identify contiguous candidate regions
            valley_mask = [val < valley_threshold for val in smoothed]
            candidate_regions: List[CandidateRegion] = []
            recommended_positions: List[int] = []

            in_region = False
            r_start = 0

            for idx, is_v in enumerate(valley_mask):
                pos = idx + 1  # 1-based
                if is_v:
                    recommended_positions.append(pos)
                    if not in_region:
                        in_region = True
                        r_start = pos
                else:
                    if in_region:
                        in_region = False
                        r_end = pos - 1
                        r_len = r_end - r_start + 1
                        if r_len >= min_region_len:
                            candidate_regions.append(
                                CandidateRegion(
                                    start_1based=r_start,
                                    end_1based=r_end,
                                    length=r_len,
                                    region_type="predicted_flexible_hinge",
                                    reason=f"Contact density valley (< {valley_threshold:.2f}, baseline mean={mean_c:.2f})",
                                    heuristic_basis=f"ESM2_{self.model_name}_non_local_contact_density",
                                )
                            )

            # Handle trailing region
            if in_region:
                r_end = n_res
                r_len = r_end - r_start + 1
                if r_len >= min_region_len:
                    candidate_regions.append(
                        CandidateRegion(
                            start_1based=r_start,
                            end_1based=r_end,
                            length=r_len,
                            region_type="predicted_flexible_hinge",
                            reason=f"Contact density valley (< {valley_threshold:.2f}, baseline mean={mean_c:.2f})",
                            heuristic_basis=f"ESM2_{self.model_name}_non_local_contact_density",
                        )
                    )

            coverage = round(len(recommended_positions) / float(n_res), 4)

            return ESMRecommendationReport(
                protein_length=n_res,
                model_name=self.model_name,
                contact_density_min=round(c_min, 3),
                contact_density_mean=round(mean_c, 3),
                contact_density_max=round(c_max, 3),
                candidate_regions=candidate_regions,
                recommended_positions=recommended_positions,
                coverage=coverage,
                status="AVAILABLE",
                reviewed=False,  # Hard invariant
            )

        except Exception as exc:
            return ESMRecommendationReport(
                protein_length=n_res,
                model_name=self.model_name,
                contact_density_min=0.0,
                contact_density_mean=0.0,
                contact_density_max=0.0,
                status="ERROR",
                reason=f"ESM-2 inference error: {exc}",
            )


def recommend_esm_pause_regions(
    protein_sequence: str,
    model_name: str = "esm2_t6_8M_UR50D",
    threshold_std: float = 0.75,
    min_region_len: int = 3,
) -> ESMRecommendationReport:
    """Convenience functional interface for ESM-2 structural pause recommendations."""
    annotator = ESMStructureAnnotator(model_name=model_name)
    return annotator.annotate(
        protein_sequence=protein_sequence,
        threshold_std=threshold_std,
        min_region_len=min_region_len,
    )
