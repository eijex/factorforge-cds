# factorforge/src/factorforge/rules/registry.py
"""Declarative Rule Registry for FactorForge bio-compiler.

Provides unified rule governance across biological invariants, assembly constraints,
synthesis risks, and regulatory recommendations.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Dict, List, Optional, Set

from factorforge.rules.models import (
    AuthorityType,
    EnforcementLevel,
    RuleAuthority,
    RuleCategory,
    RuleDefinition,
    RuleScope,
)


# ---------------------------------------------------------------------------
# Evaluator Implementations
# ---------------------------------------------------------------------------

def _check_bsai(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """BsaI Type IIS recognition site check (GGTCTC / GAGACC)."""
    seq_upper = sequence.upper()
    forward_matches = [m.start() for m in re.finditer(r"GGTCTC", seq_upper)]
    reverse_matches = [m.start() for m in re.finditer(r"GAGACC", seq_upper)]
    total_sites = len(forward_matches) + len(reverse_matches)
    return {
        "passed": total_sites == 0,
        "site_count": total_sites,
        "forward_positions": forward_matches,
        "reverse_positions": reverse_matches,
    }


def _check_bsmbi(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """BsmBI / Esp3I Type IIS recognition site check (CGTCTC / GAGACG)."""
    seq_upper = sequence.upper()
    forward_matches = [m.start() for m in re.finditer(r"CGTCTC", seq_upper)]
    reverse_matches = [m.start() for m in re.finditer(r"GAGACG", seq_upper)]
    total_sites = len(forward_matches) + len(reverse_matches)
    return {
        "passed": total_sites == 0,
        "site_count": total_sites,
        "forward_positions": forward_matches,
        "reverse_positions": reverse_matches,
    }


def _check_reading_frame(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """CDS length must be an exact multiple of 3 and start with ATG."""
    seq_upper = sequence.upper()
    is_multiple_3 = len(seq_upper) % 3 == 0
    starts_with_atg = seq_upper.startswith("ATG") if len(seq_upper) >= 3 else False
    passed = is_multiple_3 and starts_with_atg
    return {
        "passed": passed,
        "length": len(seq_upper),
        "is_multiple_of_three": is_multiple_3,
        "starts_with_atg": starts_with_atg,
    }


def _check_cryptic_splice(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """Plant/Eukaryotic canonical splice donor consensus motif ([A/C]AGGT[A/G]AGT)."""
    seq_upper = sequence.upper()
    donor_matches = [m.start() for m in re.finditer(r"[AC]AGGT[AG]AGT", seq_upper)]
    return {
        "passed": len(donor_matches) == 0,
        "risk_count": len(donor_matches),
        "donor_positions": donor_matches,
    }


def _check_polya_signals(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """Premature polyadenylation signals (AATAAA, ATTAAA, AATAAT)."""
    seq_upper = sequence.upper()
    polya_matches = [m.start() for m in re.finditer(r"(AATAAA|ATTAAA|AATAAT)", seq_upper)]
    return {
        "passed": len(polya_matches) == 0,
        "count": len(polya_matches),
        "positions": polya_matches,
    }


def _check_au_rich_elements(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """AU-rich mRNA destabilizing elements (ATTTA / AUUUA)."""
    seq_upper = sequence.upper()
    are_matches = [m.start() for m in re.finditer(r"ATTTA", seq_upper)]
    return {
        "passed": len(are_matches) == 0,
        "count": len(are_matches),
        "positions": are_matches,
    }


def _check_homopolymer_runs(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """Flag homopolymer runs >= 6 nt (synthesis yield risk)."""
    seq_upper = sequence.upper()
    runs = [m.span() for m in re.finditer(r"(A{6,}|C{6,}|G{6,}|T{6,})", seq_upper)]
    return {
        "passed": len(runs) == 0,
        "run_count": len(runs),
        "run_spans": runs,
    }


def _check_gc_extremes(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """Sliding window local GC extremes (window=50, threshold: <25% or >75%)."""
    seq_upper = sequence.upper()
    w = 50
    extreme_windows = []
    if len(seq_upper) >= w:
        for i in range(len(seq_upper) - w + 1):
            sub = seq_upper[i : i + w]
            gc = (sub.count("G") + sub.count("C")) / w
            if gc < 0.25 or gc > 0.75:
                extreme_windows.append((i, round(gc, 3)))
    return {
        "passed": len(extreme_windows) == 0,
        "extreme_count": len(extreme_windows),
        "extreme_windows": extreme_windows[:10],
    }


def _check_initiation_mfe_policy(sequence: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """Advisory check on 5' translation initiation window structural risk."""
    initiation_mfe = context.get("initiation_mfe")
    if initiation_mfe is None:
        return {"passed": True, "details": "initiation_mfe metric not provided in context"}
    # Advisory threshold: MFE should ideally be >= -3.0 kcal/mol
    passed = float(initiation_mfe) >= -3.0
    return {
        "passed": passed,
        "initiation_mfe": initiation_mfe,
        "threshold": -3.0,
    }


# ---------------------------------------------------------------------------
# Rule Registry Class
# ---------------------------------------------------------------------------

class RuleRegistry:
    """Registry managing collection of versioned, scoped, and attributed rules."""

    def __init__(self, ruleset_id: str = "default-assembly-2026.09"):
        self.ruleset_id = ruleset_id
        self._rules: Dict[str, RuleDefinition] = {}
        self._register_default_rules()

    def _register_default_rules(self) -> None:
        """Register the standard suite across Hard Invariants, Warnings, and Policies."""
        
        # 1. Hard Invariants (hard_fail)
        self.register(
            RuleDefinition(
                rule_id="assembly.type_iis.bsai.v1",
                name="BsaI Type IIS Absence",
                description="Strict exclusion of BsaI recognition sites (GGTCTC/GAGACC) for GoldenBraid / GoldenGate assembly.",
                category=RuleCategory.ASSEMBLY,
                enforcement=EnforcementLevel.HARD_FAIL,
                authority=RuleAuthority(
                    authority_type=AuthorityType.EIJEX_INTERNAL_POLICY,
                    source_name="assembly-profile-goldenbraid-v2",
                ),
                scope=RuleScope(
                    assembly_methods={"golden_gate", "golden_braid", "*"},
                    molecule_types={"dna"},
                ),
                version="1.0.0",
                evaluator_fn=_check_bsai,
            )
        )

        self.register(
            RuleDefinition(
                rule_id="assembly.type_iis.bsmbi.v1",
                name="BsmBI / Esp3I Type IIS Absence",
                description="Strict exclusion of BsmBI recognition sites (CGTCTC/GAGACG) for standard Level 2 GoldenGate assembly.",
                category=RuleCategory.ASSEMBLY,
                enforcement=EnforcementLevel.HARD_FAIL,
                authority=RuleAuthority(
                    authority_type=AuthorityType.EIJEX_INTERNAL_POLICY,
                    source_name="assembly-profile-goldengate-v1",
                ),
                scope=RuleScope(
                    assembly_methods={"golden_gate", "*"},
                    molecule_types={"dna"},
                ),
                version="1.0.0",
                evaluator_fn=_check_bsmbi,
            )
        )

        self.register(
            RuleDefinition(
                rule_id="biological.reading_frame.v1",
                name="CDS Reading Frame & Start Codon",
                description="CDS must have a length multiple of 3 and start with a valid ATG initiation codon.",
                category=RuleCategory.ASSEMBLY,
                enforcement=EnforcementLevel.HARD_FAIL,
                authority=RuleAuthority(
                    authority_type=AuthorityType.COMMUNITY_STANDARD,
                    source_name="standard_genetic_code",
                ),
                scope=RuleScope(molecule_types={"dna", "mrna"}),
                version="1.0.0",
                evaluator_fn=_check_reading_frame,
            )
        )

        # 2. Advisory Biological & RNA Risks (warning)
        self.register(
            RuleDefinition(
                rule_id="rna.cryptic_splice.v1",
                name="Eukaryotic Cryptic Splice Donor Risk",
                description="Advisory scan for canonical eukaryotic splice donor consensus motifs.",
                category=RuleCategory.RNA_RISK,
                enforcement=EnforcementLevel.WARNING,
                authority=RuleAuthority(
                    authority_type=AuthorityType.COMMUNITY_STANDARD,
                    source_name="plant_mrna_processing_standards",
                ),
                scope=RuleScope(
                    target_hosts={"Nicotiana benthamiana", "Arabidopsis thaliana", "Homo sapiens", "*"},
                    molecule_types={"dna", "mrna"},
                ),
                version="1.0.0",
                evaluator_fn=_check_cryptic_splice,
            )
        )

        self.register(
            RuleDefinition(
                rule_id="rna.polya_motifs.v1",
                name="Premature Polyadenylation Motif Risk",
                description="Advisory scan for premature polyadenylation consensus sequences (AATAAA/ATTAAA).",
                category=RuleCategory.RNA_RISK,
                enforcement=EnforcementLevel.WARNING,
                authority=RuleAuthority(
                    authority_type=AuthorityType.COMMUNITY_STANDARD,
                    source_name="eukaryotic_polyadenylation_guidelines",
                ),
                scope=RuleScope(molecule_types={"dna", "mrna"}),
                version="1.0.0",
                evaluator_fn=_check_polya_signals,
            )
        )

        self.register(
            RuleDefinition(
                rule_id="rna.au_rich_elements.v1",
                name="AU-Rich Destabilizing Elements (ARE)",
                description="Advisory scan for pentameric ATTTA motifs that mediate mRNA instability in eukaryotic hosts.",
                category=RuleCategory.RNA_RISK,
                enforcement=EnforcementLevel.WARNING,
                authority=RuleAuthority(
                    authority_type=AuthorityType.COMMUNITY_STANDARD,
                    source_name="mrna_decay_pathways",
                ),
                scope=RuleScope(molecule_types={"dna", "mrna"}),
                version="1.0.0",
                evaluator_fn=_check_au_rich_elements,
            )
        )

        self.register(
            RuleDefinition(
                rule_id="synthesis.gc_extremes.v1",
                name="Local GC Extremes (Sliding Window)",
                description="Flags local windows of 50 nt with GC < 25% or > 75% causing chemical synthesis dropout.",
                category=RuleCategory.INTERNAL_POLICY,
                enforcement=EnforcementLevel.WARNING,
                authority=RuleAuthority(
                    authority_type=AuthorityType.EIJEX_INTERNAL_POLICY,
                    source_name="eijex_synthesis_feasibility_guidelines",
                ),
                scope=RuleScope(molecule_types={"dna"}),
                version="1.0.0",
                evaluator_fn=_check_gc_extremes,
            )
        )

        # 3. Policy & Quality Recommendations (informational)
        self.register(
            RuleDefinition(
                rule_id="policy.synthesis.homopolymer.v1",
                name="Homopolymer Run Length Limit",
                description="Eijex policy recommendation to avoid mono-nucleotide runs >= 6 nt for gene synthesis yield.",
                category=RuleCategory.INTERNAL_POLICY,
                enforcement=EnforcementLevel.INFORMATIONAL,
                authority=RuleAuthority(
                    authority_type=AuthorityType.EIJEX_INTERNAL_POLICY,
                    source_name="eijex_gene_synthesis_standards_2026",
                ),
                scope=RuleScope(molecule_types={"dna"}),
                version="1.0.0",
                evaluator_fn=_check_homopolymer_runs,
            )
        )

        self.register(
            RuleDefinition(
                rule_id="policy.initiation_mfe.v1",
                name="5' Translation Initiation Structural Relaxation",
                description="Quality indicator for 5' initiation region secondary structure (target MFE >= -3.0 kcal/mol).",
                category=RuleCategory.INTERNAL_POLICY,
                enforcement=EnforcementLevel.INFORMATIONAL,
                authority=RuleAuthority(
                    authority_type=AuthorityType.EIJEX_INTERNAL_POLICY,
                    source_name="eijex_dp_v2_1_initiation_spec",
                ),
                scope=RuleScope(molecule_types={"dna", "mrna"}),
                version="1.0.0",
                evaluator_fn=_check_initiation_mfe_policy,
            )
        )

    def register(self, rule: RuleDefinition) -> None:
        """Register or update a rule definition."""
        self._rules[rule.rule_id] = rule

    def get_rule(self, rule_id: str) -> Optional[RuleDefinition]:
        """Retrieve rule by ID."""
        return self._rules.get(rule_id)

    def list_rules(
        self,
        assembly_method: Optional[str] = None,
        host: Optional[str] = None,
        enforcement: Optional[EnforcementLevel] = None,
    ) -> List[RuleDefinition]:
        """Filter rules by context and enforcement level."""
        results = []
        for rule in self._rules.values():
            if not rule.scope.matches(assembly_method=assembly_method, host=host):
                continue
            if enforcement and rule.enforcement != enforcement:
                continue
            results.append(rule)
        return results

    def evaluate_sequence(
        self,
        sequence: str,
        assembly_method: Optional[str] = None,
        host: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Evaluate all applicable rules and aggregate results by enforcement level."""
        active_rules = self.list_rules(assembly_method=assembly_method, host=host)
        ctx = context or {}

        hard_fails = []
        warnings = []
        info_items = []
        all_passed = True

        for rule in active_rules:
            res = rule.evaluate(sequence, ctx)
            passed = res.get("passed", True)
            entry = {
                "rule_id": rule.rule_id,
                "name": rule.name,
                "enforcement": rule.enforcement.value,
                "authority": rule.authority.authority_type.value,
                "result": res,
            }

            if not passed:
                if rule.enforcement == EnforcementLevel.HARD_FAIL:
                    hard_fails.append(entry)
                    all_passed = False
                elif rule.enforcement == EnforcementLevel.WARNING:
                    warnings.append(entry)
                else:
                    info_items.append(entry)

        return {
            "ruleset_id": self.ruleset_id,
            "ruleset_digest": self.compute_digest(),
            "all_passed": all_passed,
            "hard_fail_count": len(hard_fails),
            "warning_count": len(warnings),
            "hard_fails": hard_fails,
            "warnings": warnings,
            "informational": info_items,
        }

    def compute_digest(self) -> str:
        """Compute cryptographic SHA-256 fingerprint of the current rule registry."""
        rule_summaries = sorted([
            f"{r.rule_id}:{r.version}:{r.enforcement.value}:{r.authority.authority_type.value}"
            for r in self._rules.values()
        ])
        raw_text = ";".join(rule_summaries).encode("utf-8")
        return f"sha256:{hashlib.sha256(raw_text).hexdigest()}"
