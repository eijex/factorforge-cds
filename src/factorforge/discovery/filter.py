"""Hard Constraint Filter for FactorForge Discovery Slate (Job 283A).

Enforces non-negotiable binary biological and construct feasibility rules.
Evaluated metrics (MFE, 50-bp GC local fluctuations, repeat scores) are NOT
treated as hard fails.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set

from factorforge.analysis.metrics import (
    calculate_gc,
    detect_homopolymers,
    translate_dna,
)
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs


@dataclass
class HardConstraintConfig:
    """Configuration for binary hard constraint validation."""

    global_gc_min: float = 30.0
    global_gc_max: float = 55.0
    initiation_45nt_gc_min: float = 15.0
    initiation_45nt_gc_max: float = 45.0
    homopolymer_max_run: int = 5
    forbidden_type_iis_enzymes: Set[str] = field(
        default_factory=lambda: {"BsaI", "BsmBI", "BpiI", "SapI"}
    )
    forbidden_other_motifs: List[str] = field(
        default_factory=lambda: ["GCGGCCGC"]  # NotI
    )


@dataclass
class HardFilterResult:
    """Result of hard constraint validation on a candidate sequence."""

    is_feasible: bool
    violations: List[str]
    details: Dict[str, Any]


class HardConstraintFilter:
    """Evaluates candidate sequences against strict biological and assembly invariants."""

    def __init__(self, config: Optional[HardConstraintConfig] = None) -> None:
        self.config = config or HardConstraintConfig()

        # Build canonical motif list (Type IIS + other specific motifs)
        self.forbidden_motifs = get_canonical_forbidden_motifs(
            self.config.forbidden_type_iis_enzymes
        )
        for m in self.config.forbidden_other_motifs:
            m_upper = m.upper()
            if m_upper not in self.forbidden_motifs:
                self.forbidden_motifs.append(m_upper)

    def evaluate(
        self,
        sequence_dna: str,
        expected_protein_aa: str,
        expected_stop: Optional[str] = None,
        upstream_flank: str = "",
        downstream_flank: str = "",
    ) -> HardFilterResult:
        """Validates sequence against all binary hard constraints."""
        seq = sequence_dna.strip().upper()
        violations: List[str] = []
        details: Dict[str, Any] = {}

        # 1. Reading Frame & Length Integrity
        if len(seq) % 3 != 0:
            violations.append(f"Sequence length ({len(seq)}) is not a multiple of 3")

        # 2. Stop Codon & Exact Amino Acid Translation
        expected_clean = expected_protein_aa.strip().upper().rstrip("*")
        exp_stop = expected_stop.strip().upper() if expected_stop else None

        if exp_stop:
            if not seq.endswith(exp_stop):
                terminal_codon = seq[-3:] if len(seq) >= 3 else seq
                violations.append(
                    f"Sequence does not end with expected stop codon '{exp_stop}', found '{terminal_codon}'"
                )
            coding_dna = (
                seq[: -len(exp_stop)]
                if seq.endswith(exp_stop)
                else (seq[:-3] if len(seq) >= 3 and seq[-3:] in ("TAA", "TAG", "TGA") else seq)
            )
        else:
            if len(seq) >= 3 and seq[-3:] in ("TAA", "TAG", "TGA"):
                coding_dna = seq[:-3]
            else:
                coding_dna = seq

        translated = translate_dna(coding_dna)

        # Check internal stop codons
        if "*" in translated:
            violations.append(
                f"Internal stop codons detected in coding sequence: {translated.count('*')}"
            )

        # Check exact amino acid translation
        clean_translated = translated.rstrip("*")
        if clean_translated != expected_clean:
            violations.append(
                f"Translated amino acid sequence does not match expected target: "
                f"expected len {len(expected_clean)}, got len {len(clean_translated)}"
            )

        # 3. Full Construct Flank Junction & Forbidden Restriction Motifs (Type IIS + NotI)
        full_construct = upstream_flank.strip().upper() + seq + downstream_flank.strip().upper()
        detected_motifs: List[str] = []
        for motif in self.forbidden_motifs:
            if motif in full_construct:
                detected_motifs.append(motif)
        if detected_motifs:
            violations.append(
                f"Forbidden restriction motifs found in construct: {', '.join(detected_motifs)}"
            )
        details["forbidden_motifs_detected"] = detected_motifs

        # 4. Global GC Range Contract
        gc_percent = calculate_gc(seq)
        details["global_gc_percent"] = gc_percent
        if not (self.config.global_gc_min - 1e-4 <= gc_percent <= self.config.global_gc_max + 1e-4):
            violations.append(
                f"Global GC {gc_percent:.1f}% outside allowed contract "
                f"[{self.config.global_gc_min}%, {self.config.global_gc_max}%]"
            )

        # 5. 5' Initiation Region GC Band Contract (handles short sequences by scaling window)
        check_len = min(45, len(seq))
        if check_len > 0:
            gc_5p = calculate_gc(seq[:check_len])
            details["gc_5p_initiation"] = gc_5p
            details["gc_5p_window_nt"] = check_len
            if not (
                self.config.initiation_45nt_gc_min <= gc_5p <= self.config.initiation_45nt_gc_max
            ):
                violations.append(
                    f"5' initiation ({check_len}-nt) GC {gc_5p:.1f}% outside allowed contract "
                    f"[{self.config.initiation_45nt_gc_min}%, {self.config.initiation_45nt_gc_max}%]"
                )

        # 6. Homopolymer Bound on Full Construct (flag runs exceeding max_run)
        homopolymers = detect_homopolymers(
            full_construct, max_run=self.config.homopolymer_max_run + 1
        )
        details["homopolymers"] = homopolymers
        if homopolymers:
            violations.append(
                f"Homopolymer runs exceeding {self.config.homopolymer_max_run}-mer: "
                f"{len(homopolymers)} occurrences"
            )

        return HardFilterResult(
            is_feasible=(len(violations) == 0),
            violations=violations,
            details=details,
        )
