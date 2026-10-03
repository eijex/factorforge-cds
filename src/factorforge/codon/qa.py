"""In-Silico Quality Assurance (QA) verification suite for generated CDS candidates."""

from __future__ import annotations

import math
from typing import Dict, List, Optional
from pydantic import BaseModel
from factorforge.codon.policy import CODON_TO_AA

# Type IIS Restriction Enzyme Recognition Sequences
ENZYME_RECOGNITION = {
    "BsaI": ["GGTCTC", "GAGACC"],
    "BpiI": ["GAAGAC", "GTCTTC"],
    "BsmBI": ["CGTCTC", "GAGACG"],
}


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(complement)[::-1]


def translate(seq: str) -> str:
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    return "".join(CODON_TO_AA.get(c.upper(), "?") for c in codons)


class InSilicoQAReport(BaseModel):
    """Detailed QA evaluation report for a CDS construct."""

    construct_id: str
    length_nt: int
    length_aa: int
    translation_match: bool
    internal_stop_count: int
    terminal_stop_present: bool
    signal_peptide_preserved: bool
    type_iis_findings: Dict[str, List[int]]  # Enzyme name -> 1-based start positions
    type_iis_clean: bool
    gc_percent: float
    cai: Optional[float] = None
    cai_status: str = "unavailable"
    cai_reason: Optional[str] = None
    consecutive_repeat_count: int  # count of adjacent identical codons for same AA
    overall_passed: bool
    rejection_reasons: List[str]


class InSilicoQA:
    """Evaluates candidate sequences against computational validity rules."""

    @staticmethod
    def evaluate(
        cds: str,
        expected_protein: str,
        construct_id: str = "candidate",
        baseline_sp_cds: Optional[str] = None,
        codon_weights: Optional[Dict[str, float]] = None,
    ) -> InSilicoQAReport:
        cds = cds.upper()
        expected_protein = expected_protein.upper()

        rejection_reasons: List[str] = []

        # 1. Translation match
        actual_prot = translate(cds)
        clean_actual = actual_prot.rstrip("*")
        clean_expected = expected_protein.rstrip("*")
        translation_match = clean_actual == clean_expected
        if not translation_match:
            rejection_reasons.append(
                f"Translation mismatch: expected {len(clean_expected)} aa, got {len(clean_actual)} aa"
            )

        # 2. Stop codons
        codons = [cds[i:i+3] for i in range(0, len(cds) - 2, 3)]
        internal_stops = [i for i, c in enumerate(codons[:-1]) if CODON_TO_AA.get(c) == "*"]
        internal_stop_count = len(internal_stops)
        if internal_stop_count > 0:
            rejection_reasons.append(f"Contains {internal_stop_count} internal stop codons at indices {internal_stops}")

        terminal_stop_present = len(codons) > 0 and CODON_TO_AA.get(codons[-1]) == "*"

        # 3. Signal peptide check
        signal_peptide_preserved = True
        if baseline_sp_cds:
            sp_len = len(baseline_sp_cds)
            if cds[:sp_len] != baseline_sp_cds.upper():
                signal_peptide_preserved = False
                rejection_reasons.append(f"5' Signal peptide differs from baseline {sp_len} nt")

        # 4. Type IIS restriction sites
        type_iis_findings: Dict[str, List[int]] = {}
        type_iis_clean = True
        for enzyme, patterns in ENZYME_RECOGNITION.items():
            matches: List[int] = []
            for pat in patterns:
                start = 0
                while True:
                    idx = cds.find(pat, start)
                    if idx == -1:
                        break
                    matches.append(idx + 1)  # 1-based coordinate
                    start = idx + 1
            if matches:
                type_iis_findings[enzyme] = sorted(matches)
                type_iis_clean = False
                rejection_reasons.append(f"Found {len(matches)} {enzyme} site(s) at {matches}")
            else:
                type_iis_findings[enzyme] = []

        # 5. GC content
        gc_count = cds.count("G") + cds.count("C")
        gc_percent = round((gc_count / len(cds)) * 100, 2) if len(cds) > 0 else 0.0

        # 6. CAI calculation
        if codon_weights and len(codons) > 0:
            log_sum = sum(math.log(max(codon_weights.get(c, 0.01), 1e-4)) for c in codons)
            cai = round(math.exp(log_sum / len(codons)), 4)
            cai_status = "calculated"
            cai_reason = None
        else:
            cai = None
            cai_status = "unavailable"
            cai_reason = "host_codon_weight_reference_not_provided"

        # 7. Consecutive codon repeats
        consecutive_repeats = 0
        for i in range(len(codons) - 1):
            if codons[i] == codons[i+1]:
                consecutive_repeats += 1

        overall_passed = (
            translation_match
            and internal_stop_count == 0
            and signal_peptide_preserved
            and type_iis_clean
        )

        return InSilicoQAReport(
            construct_id=construct_id,
            length_nt=len(cds),
            length_aa=len(clean_actual),
            translation_match=translation_match,
            internal_stop_count=internal_stop_count,
            terminal_stop_present=terminal_stop_present,
            signal_peptide_preserved=signal_peptide_preserved,
            type_iis_findings=type_iis_findings,
            type_iis_clean=type_iis_clean,
            gc_percent=gc_percent,
            cai=cai,
            cai_status=cai_status,
            cai_reason=cai_reason,
            consecutive_repeat_count=consecutive_repeats,
            overall_passed=overall_passed,
            rejection_reasons=rejection_reasons,
        )
