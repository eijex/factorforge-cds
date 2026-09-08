"""Experimental FactorForge-LM decoding constraints.

This module is an inactive research scaffold. It provides deterministic token-level
masking helpers that can be tested without model weights, but it is not a public
benchmark engine and it does not make biological, synthesis, cloning, or wet-lab
performance guarantees.
"""

from __future__ import annotations

from typing import Optional

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.engines.lm.tokenizer.control_tokenizer import (
    ALL_CODONS,
    FactorForgeControlTokenizer,
)

TYPE_IIS_PATTERNS = (
    "GGTCTC",  # BsaI
    "GAGACC",  # BsaI reverse complement
    "GAAGAC",  # BpiI/BbsI
    "GTCTTC",  # BpiI/BbsI reverse complement
    "CGTCTC",  # BsmBI/Esp3I
    "GAGACG",  # BsmBI/Esp3I reverse complement
)


class FactorForgeLogitMasker:
    """Apply deterministic next-codon masks for experimental LM decoding.

    The masker is intentionally narrow: it can prevent candidate codons that would
    create configured Type IIS recognition sites and, when the expected amino acid
    is supplied, masks non-synonymous codons. GC handling is advisory/penalty-only
    and must not be interpreted as a hard design acceptance predicate.
    """

    def __init__(
        self,
        tokenizer: FactorForgeControlTokenizer,
        enable_type_iis_masking: bool = True,
        target_gc_min: float = 0.40,
        target_gc_max: float = 0.47,
    ) -> None:
        self.tokenizer = tokenizer
        self.enable_type_iis_masking = enable_type_iis_masking
        self.target_gc_min = target_gc_min
        self.target_gc_max = target_gc_max
        self.codon_to_id = {
            codon: tokenizer.token_to_id[codon]
            for codon in ALL_CODONS
            if codon in tokenizer.token_to_id
        }

    def apply_logit_masks_with_reasons(
        self,
        current_cds_prefix: str,
        logits: list[float],
        expected_next_aa: Optional[str] = None,
    ) -> tuple[list[float], dict[str, str]]:
        if len(logits) < self.tokenizer.vocab_size:
            raise ValueError("logits length is smaller than tokenizer vocabulary")

        prefix = "".join(current_cds_prefix.upper().replace("U", "T").split())
        expected = expected_next_aa.upper() if expected_next_aa else None
        masked_logits = list(logits)
        rejections = {}

        for codon, token_id in self.codon_to_id.items():
            if expected is not None and STANDARD_GENETIC_CODE.get(codon) != expected:
                masked_logits[token_id] = -float("inf")
                rejections[codon] = "HARD_REJECT: non-synonymous"
                continue

            candidate_prefix = prefix + codon
            if self.enable_type_iis_masking and any(
                pattern in candidate_prefix for pattern in TYPE_IIS_PATTERNS
            ):
                masked_logits[token_id] = -float("inf")
                rejections[codon] = "HARD_REJECT: Type IIS"
                continue

            if len(candidate_prefix) >= 30:
                gc_ratio = (candidate_prefix.count("G") + candidate_prefix.count("C")) / len(
                    candidate_prefix
                )
                codon_gc = (codon.count("G") + codon.count("C")) / 3.0
                if gc_ratio < self.target_gc_min and codon_gc == 0:
                    masked_logits[token_id] -= 5.0
                    rejections[codon] = f"SOFT_PENALTY: cumulative GC lower bound (cumulative GC {gc_ratio*100:.1f}%)"
                elif gc_ratio > self.target_gc_max and codon_gc > 0.66:
                    masked_logits[token_id] -= 5.0
                    rejections[codon] = f"SOFT_PENALTY: cumulative GC upper bound (cumulative GC {gc_ratio*100:.1f}%)"

        return masked_logits, rejections

    def apply_logit_masks(
        self,
        current_cds_prefix: str,
        logits: list[float],
        expected_next_aa: Optional[str] = None,
    ) -> list[float]:
        return self.apply_logit_masks_with_reasons(current_cds_prefix, logits, expected_next_aa)[0]
