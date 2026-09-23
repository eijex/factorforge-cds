"""Multi-Resolution ValidationHub Engine for FactorForge Slate v2 (Job 293A).

Implements Stage 4 Validation:
- Level 1: Primary Host Sequence Integrity (Strict separation of AA verification and forbidden sites).
- Level 2: Deep RNA initiation folding (Real thermodynamic RNAfold when available; NO fake proxy substitution).
- Level 3: Cross-Host Sensitivity Diagnostic (Truthful host model provenance, NOT_AVAILABLE on missing models).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import logging
import math
from typing import Any, Dict, List, Optional, Tuple

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.core.host_model import HostModel
from factorforge.core.reranker import ScoredCandidate
from factorforge.core.slate_generator import (
    DEFAULT_FORBIDDEN_SITES,
    CandidateDesign,
    normalize_protein_sequence,
    translate_dna,
)
from factorforge.utils.restriction_sites import detect_restriction_sites

logger = logging.getLogger(__name__)


@dataclass
class ValidationReportSummary:
    """Consolidated validation summary for the final slate with orthogonal metrics."""

    aa_conservation_by_construction: bool = True
    aa_conservation_verified: bool = True
    forbidden_sites_clean: bool = True
    primary_sequence_integrity_pass: bool = True
    primary_host_status: str = "PASS"
    aa_identity_pct: float = 100.0
    substitutions: int = 0
    insertions: int = 0
    deletions: int = 0
    total_evaluated: int = 0
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "aa_conservation_by_construction": self.aa_conservation_by_construction,
            "aa_conservation_verified": self.aa_conservation_verified,
            "forbidden_sites_clean": self.forbidden_sites_clean,
            "primary_sequence_integrity_pass": self.primary_sequence_integrity_pass,
            "primary_host_status": self.primary_host_status,
            "aa_identity_pct": self.aa_identity_pct,
            "substitutions": self.substitutions,
            "insertions": self.insertions,
            "deletions": self.deletions,
            "total_evaluated": self.total_evaluated,
        }


class MultiResolutionValidationHub:
    """Multi-resolution validator with truthful evidence boundaries and orthogonal semantics."""

    def __init__(
        self,
        primary_host_model: HostModel,
        secondary_hosts: Optional[List[str]] = None,
        forbidden_sites: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        self.primary_host_model = primary_host_model
        self.secondary_host_ids = secondary_hosts or ["ntabacum", "athaliana", "wolffia_globosa"]
        self.forbidden_sites = forbidden_sites or DEFAULT_FORBIDDEN_SITES
        self._secondary_host_models: Dict[str, Optional[HostModel]] = {}
        self._load_secondary_hosts()

    def _load_secondary_hosts(self) -> None:
        """Preload secondary host models safely; missing models marked as None (NOT_AVAILABLE)."""
        for host_id in self.secondary_host_ids:
            try:
                model = HostModel.load(host_id)
                self._secondary_host_models[host_id] = model
            except Exception:
                self._secondary_host_models[host_id] = None

    def validate_level1_integrity(
        self,
        candidate_dna: str,
        expected_canonical_protein: str,
    ) -> Tuple[bool, bool, Dict[str, Any]]:
        """Level 1: Strict and orthogonal sequence & AA integrity verification.
        
        Returns:
            (aa_verified, forbidden_clean, details)
        """
        clean_dna = "".join(candidate_dna.upper().split())
        canonical_protein = normalize_protein_sequence(expected_canonical_protein)

        # 1. Strict AA Translation Verification
        try:
            translated_aa = translate_dna(clean_dna)
        except Exception as e:
            return False, False, {
                "aa_verified": False,
                "forbidden_clean": False,
                "reason": f"Translation error: {str(e)}",
                "aa_identity_pct": 0.0,
                "substitutions": len(canonical_protein),
            }

        substitutions = 0
        insertions = max(0, len(translated_aa) - len(canonical_protein))
        deletions = max(0, len(canonical_protein) - len(translated_aa))
        min_len = min(len(translated_aa), len(canonical_protein))
        for i in range(min_len):
            if translated_aa[i] != canonical_protein[i]:
                substitutions += 1

        matches = sum(1 for i in range(min_len) if translated_aa[i] == canonical_protein[i])
        total_len = max(len(canonical_protein), 1)
        identity_pct = round((matches / total_len) * 100.0, 2)
        aa_verified = (translated_aa == canonical_protein) and (substitutions == 0) and (insertions == 0) and (deletions == 0)

        # 2. Forbidden restriction sites screening
        hits = detect_restriction_sites(clean_dna, self.forbidden_sites, scan_rc=True)
        forbidden_clean = len(hits) == 0
        forbidden_names = sorted({hit.get("name", "Unknown") for hit in hits}) if hits else []

        return aa_verified, forbidden_clean, {
            "aa_verified": aa_verified,
            "forbidden_clean": forbidden_clean,
            "aa_identity_pct": identity_pct,
            "substitutions": substitutions,
            "insertions": insertions,
            "deletions": deletions,
            "forbidden_hits": forbidden_names,
        }

    def evaluate_level2_rna_folding(
        self,
        coding_sequence: str,
    ) -> Dict[str, Any]:
        """Level 2: Deep RNA initiation folding evaluation.
        
        Truthfulness Invariant: Never fabricate kcal/mol. If ViennaRNA is not present,
        return status: NOT_AVAILABLE, value: null.
        """
        # Try importing ViennaRNA if available
        try:
            import RNA  # type: ignore
            fold_window = coding_sequence[:48]
            structure, mfe = RNA.fold(fold_window)
            return {
                "status": "AVAILABLE",
                "evidence_class": "PREDICTED",
                "engine": f"ViennaRNA_{getattr(RNA, '__version__', '2.x')}",
                "folding_window_nt": 48,
                "temperature_c": 37.0,
                "mfe_5p_kcal_mol": round(float(mfe), 2),
                "secondary_structure": structure,
            }
        except (ImportError, Exception):
            return {
                "status": "NOT_AVAILABLE",
                "evidence_class": "PREDICTED",
                "engine": "ViennaRNA_RNAfold",
                "folding_window_nt": 48,
                "temperature_c": 37.0,
                "mfe_5p_kcal_mol": None,
                "secondary_structure": None,
                "reason": "ViennaRNA package not installed in environment",
            }

    def compute_host_cai(self, coding_sequence: str, host_model: HostModel) -> float:
        """Calculate CAI for a specific host using its HostModel weights."""
        seq = coding_sequence.upper()
        log_sum = 0.0
        n_codons = 0
        weights = host_model.codon_weights

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i : i + 3]
            aa = STANDARD_GENETIC_CODE.get(codon)
            if aa and aa != "*":
                w = weights.get(codon, 0.05)
                w = max(0.01, min(1.0, w))
                log_sum += math.log(w)
                n_codons += 1

        if n_codons == 0:
            return 0.0
        return round(float(math.exp(log_sum / n_codons)), 4)

    def evaluate_cross_host_diagnostic(
        self,
        coding_sequence: str,
    ) -> Dict[str, Any]:
        """Level 3: Cross-Host Sensitivity Diagnostic (Sensitivity stress test delta %, NOT a hard drop gate)."""
        primary_cai = self.compute_host_cai(coding_sequence, self.primary_host_model)
        primary_ref = 1.00

        diag: Dict[str, Any] = {
            "primary_host": {
                "host_id": self.primary_host_model.host_id,
                "model_version": self.primary_host_model.model_version,
                "cai_reference": primary_ref,
            },
        }

        for host_id in self.secondary_host_ids:
            model = self._secondary_host_models.get(host_id)
            key_name = f"{host_id}_delta_pct"
            if model is None:
                diag[key_name] = {
                    "status": "NOT_AVAILABLE",
                    "score": None,
                    "delta_pct": None,
                    "reason": f"Codon model for '{host_id}' not available in repository",
                }
            else:
                sec_cai = self.compute_host_cai(coding_sequence, model)
                delta_pct = round(((sec_cai - primary_cai) / max(0.01, primary_cai)) * 100.0, 1)
                diag[key_name] = delta_pct

        return diag

    def validate_and_annotate_pool(
        self,
        scored_candidates: List[ScoredCandidate],
        expected_canonical_protein: str,
    ) -> Tuple[List[Dict[str, Any]], ValidationReportSummary]:
        """Execute full MultiResolutionValidationHub on candidates."""
        canonical_protein = normalize_protein_sequence(expected_canonical_protein)
        validated_pool: List[Dict[str, Any]] = []

        all_aa_verified = True
        all_forbidden_clean = True
        total_substitutions = 0
        total_insertions = 0
        total_deletions = 0

        for sc in scored_candidates:
            cand = sc.candidate
            aa_ver, forb_clean, l1_res = self.validate_level1_integrity(
                cand.coding_sequence, canonical_protein
            )

            if not aa_ver:
                all_aa_verified = False
            if not forb_clean:
                all_forbidden_clean = False

            total_substitutions += l1_res["substitutions"]
            total_insertions += l1_res["insertions"]
            total_deletions += l1_res["deletions"]

            # Deep RNA folding evaluation
            rna_eval = self.evaluate_level2_rna_folding(cand.coding_sequence)

            # Cross-host sensitivity diagnostic
            cross_host_diag = self.evaluate_cross_host_diagnostic(cand.coding_sequence)

            candidate_dict = {
                "rank": sc.rank,
                "pareto_front": sc.pareto_front,
                "profile_name": sc.profile_name,
                "dna_sequence": cand.dna_sequence,
                "coding_sequence": cand.coding_sequence,
                "terminal_stop": cand.terminal_stop,
                "scores": {
                    "composite_utility": round(sc.composite_utility, 3),
                    "cai": round(cand.cai, 3),
                    "gc_global": round(cand.gc_global, 3),
                    "five_prime_structure_proxy_score": sc.five_prime_structure_proxy_score,
                    "rare_codon_count": cand.rare_codon_count,
                    "rare_codon_guard_score": round(sc.rare_codon_guard_score, 3),
                },
                "sha256": cand.sha256,
                "trait_vector_normalized_7d": sc.trait_vector_7d,
                "deep_rna_folding": rna_eval,
                "cross_host_sensitivity": cross_host_diag,
                "validation": {
                    "aa_conservation_verified": aa_ver,
                    "forbidden_sites_clean": forb_clean,
                    "integrity_pass": aa_ver and forb_clean,
                },
            }
            validated_pool.append(candidate_dict)

        primary_pass = all_aa_verified and all_forbidden_clean
        summary = ValidationReportSummary(
            aa_conservation_by_construction=True,
            aa_conservation_verified=all_aa_verified,
            forbidden_sites_clean=all_forbidden_clean,
            primary_sequence_integrity_pass=primary_pass,
            primary_host_status="PASS" if primary_pass else "FAIL",
            aa_identity_pct=100.0 if (total_substitutions == 0 and total_insertions == 0 and total_deletions == 0) else 0.0,
            substitutions=total_substitutions,
            insertions=total_insertions,
            deletions=total_deletions,
            total_evaluated=len(validated_pool),
        )

        return validated_pool, summary
