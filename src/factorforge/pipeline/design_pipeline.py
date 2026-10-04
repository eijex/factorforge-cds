"""
FactorForge Deterministic Design Pipeline (Job 355).

A pure, reproducible, non-LLM execution pipeline that transforms protein input,
frozen boundaries, and structured CodonPolicy objects into a structured DesignExperiment
containing arbitrary variants, multi-level EvidenceRecords, and deterministic auto-repair.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import re
from typing import Any, Optional

from factorforge.core.experiment import (
    ComparisonMatrix,
    DesignExperiment,
    EvidenceRecord,
    Severity,
    ValidationLevel,
    Variant,
)
from factorforge.core.policy import (
    CodonPolicy,
    PLANTFORM_BALANCED_CODON_V1,
    get_policy,
)

# Standard Genetic Code
GENETIC_CODE: dict[str, str] = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}

# Type IIS Restriction Enzymes
TYPE_IIS_PATTERNS = {
    "assembly.type_iis.bsai": (re.compile(r"GGTCTC", re.IGNORECASE), re.compile(r"GAGACC", re.IGNORECASE)),
    "assembly.type_iis.bpii": (re.compile(r"GAAGAC", re.IGNORECASE), re.compile(r"GTCTTC", re.IGNORECASE)),
    "assembly.type_iis.bsmbi": (re.compile(r"CGTCTC", re.IGNORECASE), re.compile(r"GAGACG", re.IGNORECASE)),
}


def translate_cds(dna: str) -> str:
    """Translate DNA to amino acid sequence strictly in frame."""
    dna = dna.upper().replace("U", "T")
    aa_list = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa = GENETIC_CODE.get(codon, "X")
        if aa == "*":
            break
        aa_list.append(aa)
    return "".join(aa_list)


def calculate_gc(dna: str) -> float:
    """Calculate GC percentage of a DNA sequence."""
    dna = dna.upper()
    if not dna:
        return 0.0
    gc_count = dna.count("G") + dna.count("C")
    return round((gc_count / len(dna)) * 100.0, 2)


def hamming_distance(seq1: str, seq2: str) -> int:
    """Calculate nucleotide divergence count between two equal-length sequences."""
    return sum(1 for a, b in zip(seq1, seq2) if a != b) + abs(len(seq1) - len(seq2))


class DesignPipeline:
    """
    Deterministic Design Execution Core.
    
    Guarantees:
    - 100% deterministic (no LLM, identical inputs produce bit-for-bit identical results).
    - Preserves frozen boundary (e.g. 1-60 nt) with 0% mutation.
    - Uses structured CodonPolicy object (Hamilton Largest Remainder + repeat avoidance).
    - Deterministic synonymous auto-repair for Type IIS restriction sites (BsaI, BpiI, BsmBI).
    - Motif-local, coordinate-aware cryptic splice donor remediation.
    - Emits structured DesignExperiment with arbitrary Variant instances and EvidenceRecords.
    """

    def __init__(
        self,
        host: str = "nbenthamiana",
        policy: str | CodonPolicy = "plantform_balanced_codon_v1",
        tool_version: str = "v3.5.4-job355",
    ) -> None:
        self.host = host
        self.policy: CodonPolicy = get_policy(policy)
        self.policy_name: str = self.policy.name
        self.tool_version = tool_version

    def run(
        self,
        target_name: str,
        protein_sequence: str,
        seed_dna: Optional[str] = None,
        experiment_id: str = "EXP-DESIGN-01",
        hypothesis: str = "Multi-variant experimental design",
        frozen_boundary_nt: int = 60,
        interventions: Optional[list[str]] = None,
        custom_metadata: Optional[dict[str, Any]] = None,
        fixed_timestamp: Optional[str] = None,
    ) -> DesignExperiment:
        """
        Execute the deterministic design pipeline.
        
        Args:
            target_name: Target molecule identifier (e.g. "Adalimumab Light Chain")
            protein_sequence: Full target amino acid sequence
            seed_dna: Optional seed DNA sequence (used to freeze 5' signal peptide context)
            experiment_id: Unique experiment ID
            hypothesis: Scientific intent of the experiment
            frozen_boundary_nt: Number of nucleotides at 5' end to freeze immutably (default: 60)
            interventions: List of variants to generate (default: ["baseline", "watermark"])
            custom_metadata: Extra metadata to attach
            fixed_timestamp: Optional fixed timestamp string for bit-for-bit determinism
        """
        if interventions is None:
            interventions = ["baseline", "watermark"]

        timestamp = fixed_timestamp or datetime.now(timezone.utc).isoformat()
        clean_protein = protein_sequence.strip().upper()
        frozen_regions = [
            {"start": 1, "end": frozen_boundary_nt, "unit": "nt", "annotation": "5' Signal Peptide Working Boundary"}
        ]

        experiment = DesignExperiment(
            experiment_id=experiment_id,
            target_name=target_name,
            host_organism=self.host,
            hypothesis=hypothesis,
            frozen_regions=frozen_regions,
            metadata=custom_metadata or {},
        )

        # 1. Determine frozen 5' context
        frozen_aa_len = frozen_boundary_nt // 3
        if seed_dna and len(seed_dna) >= frozen_boundary_nt:
            frozen_dna = seed_dna[:frozen_boundary_nt].upper()
            frozen_aa = translate_cds(frozen_dna)
            if frozen_aa != clean_protein[:frozen_aa_len]:
                raise ValueError(
                    f"Seed DNA 5' translation ('{frozen_aa}') does not match target protein ('{clean_protein[:frozen_aa_len]}')."
                )
        else:
            frozen_dna = self._generate_naive_sp_dna(clean_protein[:frozen_aa_len])

        mature_protein = clean_protein[frozen_aa_len:]

        # 2. Generate Baseline Variant using structured CodonPolicy (Hamilton allocation)
        baseline_mature_dna = self._allocate_mature_cds(mature_protein, watermark=False)
        baseline_raw_dna = frozen_dna + baseline_mature_dna + "TAA"
        # Deterministically repair Type IIS sites outside frozen boundary
        baseline_dna = self._repair_type_iis_sites(baseline_raw_dna, frozen_boundary_nt)
        baseline_id = f"{target_name.replace(' ', '_')}-VAR-01-BASELINE"

        baseline_variant = self._create_variant(
            variant_id=baseline_id,
            parent_variant_id=None,
            design_intent="Clean biological baseline without modification",
            interventions=["none"],
            sequence=baseline_dna,
            protein_sequence=clean_protein,
            frozen_boundary_nt=frozen_boundary_nt,
            seed_dna=seed_dna,
            timestamp=timestamp,
        )
        experiment.add_variant(baseline_variant)

        # 3. Generate Intervened Variants
        var_idx = 2
        for itype in interventions:
            if itype.lower() in ("baseline", "none"):
                continue

            if itype.lower() == "watermark":
                watermarked_mature_dna = self._allocate_mature_cds(mature_protein, watermark=True)
                raw_wm_dna = frozen_dna + watermarked_mature_dna + "TAA"
                var_dna = self._repair_type_iis_sites(raw_wm_dna, frozen_boundary_nt)
                var_id = f"{target_name.replace(' ', '_')}-VAR-{var_idx:02d}-WATERMARKED"

                var = self._create_variant(
                    variant_id=var_id,
                    parent_variant_id=baseline_id,
                    design_intent="Embedded digital watermark for sequence lineage provenance",
                    interventions=["watermark"],
                    sequence=var_dna,
                    protein_sequence=clean_protein,
                    frozen_boundary_nt=frozen_boundary_nt,
                    seed_dna=seed_dna,
                    timestamp=timestamp,
                )
                experiment.add_variant(var)

                # Pairwise comparison against baseline
                gc_shift = round(var.gc_percent - baseline_variant.gc_percent, 2)
                comp_meta: dict[str, Any] = {}
                if abs(gc_shift) > 5.0:
                    comp_meta["confounder_warning"] = (
                        f"GC shift of {gc_shift:+.2f}%p exceeds 5.0%p threshold; "
                        f"expression differential risks being confounded by mRNA secondary structure or GC content."
                    )
                comp = ComparisonMatrix(
                    comparison_name="Watermark Provenance Confounder Check",
                    variant_a_id=baseline_id,
                    variant_b_id=var_id,
                    variable_evaluated="digital watermark provenance tagging",
                    gc_shift_percent=gc_shift,
                    sequence_divergence_nt=hamming_distance(baseline_dna, var_dna),
                    metadata=comp_meta,
                )
                experiment.add_comparison(comp)
                var_idx += 1

            elif itype.lower() == "splice_remediation":
                # Coordinate-aware, motif-local synonymous remediation on baseline DNA
                remediated_dna = self._remediate_splice_sites_locally(baseline_dna, frozen_boundary_nt)
                # Ensure no Type IIS site was introduced during remediation
                remediated_dna = self._repair_type_iis_sites(remediated_dna, frozen_boundary_nt)
                var_id = f"{target_name.replace(' ', '_')}-VAR-{var_idx:02d}-SPLICE-REMEDIATED"

                var = self._create_variant(
                    variant_id=var_id,
                    parent_variant_id=baseline_id,
                    design_intent="Motif-local synonymous removal of predicted sense cryptic splice donor motifs",
                    interventions=["splice_remediation"],
                    sequence=remediated_dna,
                    protein_sequence=clean_protein,
                    frozen_boundary_nt=frozen_boundary_nt,
                    seed_dna=seed_dna,
                    timestamp=timestamp,
                )
                experiment.add_variant(var)

                comp = ComparisonMatrix(
                    comparison_name="Splice Site Remediation Impact",
                    variant_a_id=baseline_id,
                    variant_b_id=var_id,
                    variable_evaluated="synonymous cryptic splice removal",
                    gc_shift_percent=round(var.gc_percent - baseline_variant.gc_percent, 2),
                    sequence_divergence_nt=hamming_distance(baseline_dna, remediated_dna),
                )
                experiment.add_comparison(comp)
                var_idx += 1

            elif itype.lower() == "skewed_contrast":
                # Max CAI single codon bias for empirical contrast
                from factorforge.core.policy import MAX_CAI_EXTREME_V1
                skewed_mature = "".join(MAX_CAI_EXTREME_V1.allocate_sequence(mature_protein))
                raw_skewed = frozen_dna + skewed_mature + "TAA"
                var_dna = self._repair_type_iis_sites(raw_skewed, frozen_boundary_nt)
                var_id = f"{target_name.replace(' ', '_')}-VAR-{var_idx:02d}-SKEWED-CONTRAST"

                var = self._create_variant(
                    variant_id=var_id,
                    parent_variant_id=baseline_id,
                    design_intent="High-frequency biased synonymous codon allocation for empirical contrast",
                    interventions=["skewed_contrast"],
                    sequence=var_dna,
                    protein_sequence=clean_protein,
                    frozen_boundary_nt=frozen_boundary_nt,
                    seed_dna=seed_dna,
                    timestamp=timestamp,
                )
                experiment.add_variant(var)

                comp = ComparisonMatrix(
                    comparison_name="Codon Bias Skew Contrast",
                    variant_a_id=baseline_id,
                    variant_b_id=var_id,
                    variable_evaluated="codon frequency arithmetic bias",
                    gc_shift_percent=round(var.gc_percent - baseline_variant.gc_percent, 2),
                    sequence_divergence_nt=hamming_distance(baseline_dna, var_dna),
                )
                experiment.add_comparison(comp)
                var_idx += 1

        return experiment

    def _create_variant(
        self,
        variant_id: str,
        parent_variant_id: Optional[str],
        design_intent: str,
        interventions: list[str],
        sequence: str,
        protein_sequence: str,
        frozen_boundary_nt: int,
        seed_dna: Optional[str],
        timestamp: str,
    ) -> Variant:
        """Create a Variant and execute structured multi-layer evidence evaluation."""
        gc_val = calculate_gc(sequence)
        # CAI is NOT_COMPUTED until empirical / host codon reference table is linked
        cai_val: Optional[float] = None

        evidence: list[EvidenceRecord] = []

        # 1. Translation Identity (Deterministic Check)
        trans_aa = translate_cds(sequence)
        aa_match = (trans_aa == protein_sequence)
        evidence.append(
            EvidenceRecord(
                check_id="genetics.translation_identity",
                category="genetics",
                result="PASS" if aa_match else "FAIL",
                severity=Severity.HARD_FAIL,
                method="StandardGeneticCodeTranslator",
                tool_version=self.tool_version,
                confidence=1.0,
                evidence_detail=f"Translated {len(trans_aa)} aa matches expected {len(protein_sequence)} aa exactly.",
                validation_level=ValidationLevel.DETERMINISTIC_CHECK,
                timestamp=timestamp,
            )
        )

        # 2. Frozen 5' Boundary Preservation
        if seed_dna:
            expected_frozen = seed_dna[:frozen_boundary_nt].upper()
            actual_frozen = sequence[:frozen_boundary_nt].upper()
            frozen_match = (expected_frozen == actual_frozen)
            evidence.append(
                EvidenceRecord(
                    check_id="constraints.frozen_5prime_boundary",
                    category="cloning_hygiene",
                    result="PASS" if frozen_match else "FAIL",
                    severity=Severity.HARD_FAIL,
                    method="ExactPrefixComparator",
                    tool_version=self.tool_version,
                    confidence=1.0,
                    evidence_detail=f"First {frozen_boundary_nt} nt immutable boundary identical to configured seed.",
                    validation_level=ValidationLevel.DETERMINISTIC_CHECK,
                    timestamp=timestamp,
                )
            )

        # 3. Type IIS Restriction Hygiene (BsaI, BpiI, BsmBI)
        for check_id, (fwd_pat, rev_pat) in TYPE_IIS_PATTERNS.items():
            fwd_hits = len(fwd_pat.findall(sequence))
            rev_hits = len(rev_pat.findall(sequence))
            total_hits = fwd_hits + rev_hits
            evidence.append(
                EvidenceRecord(
                    check_id=check_id,
                    category="cloning_hygiene",
                    result="PASS" if total_hits == 0 else "FAIL",
                    severity=Severity.HARD_FAIL,
                    method="TypeIISRegexScanner",
                    tool_version=self.tool_version,
                    confidence=1.0,
                    evidence_detail=f"{total_hits} recognition sites detected (0 required for Golden Gate assembly).",
                    validation_level=ValidationLevel.DETERMINISTIC_CHECK,
                    timestamp=timestamp,
                )
            )

        # 4. Cryptic Splice Advisory Scan (Computational Prediction)
        splice_hits = len(re.findall(r"AGGT[AG]", sequence, re.IGNORECASE))
        if splice_hits > 0:
            evidence.append(
                EvidenceRecord(
                    check_id="rna.cryptic_splice",
                    category="rna_stability",
                    result="WARNING",
                    severity=Severity.WARNING,
                    method="SpliceDonorHeuristicScanner",
                    tool_version=self.tool_version,
                    confidence=0.72,
                    evidence_detail=f"{splice_hits} potential donor-like consensus motifs observed.",
                    validation_level=ValidationLevel.COMPUTATIONAL_PREDICTION,
                    timestamp=timestamp,
                )
            )
        else:
            evidence.append(
                EvidenceRecord(
                    check_id="rna.cryptic_splice",
                    category="rna_stability",
                    result="PASS",
                    severity=Severity.WARNING,
                    method="SpliceDonorHeuristicScanner",
                    tool_version=self.tool_version,
                    confidence=0.85,
                    evidence_detail="Zero consensus splice donor motifs detected.",
                    validation_level=ValidationLevel.COMPUTATIONAL_PREDICTION,
                    timestamp=timestamp,
                )
            )

        # 5. CAI Evaluation (Explicitly NOT_COMPUTED)
        evidence.append(
            EvidenceRecord(
                check_id="genetics.cai_score",
                category="genetics",
                result="NOT_COMPUTED",
                severity=Severity.INFO,
                method="HarmonicCodonUsageIndices",
                tool_version=self.tool_version,
                confidence=None,
                evidence_detail="Reference organism codon index table not linked; score uncomputed.",
                validation_level=ValidationLevel.NOT_VALIDATED,
                timestamp=timestamp,
            )
        )

        # 6. Wet-Lab Status (Marked NOT_VALIDATED, NOT_TESTED)
        evidence.append(
            EvidenceRecord(
                check_id="wet_lab.expression_status",
                category="wet_lab",
                result="NOT_TESTED",
                severity=Severity.INFO,
                method="InVivoPlantExpression",
                tool_version="unassigned",
                confidence=None,
                evidence_detail="Construct computationally verified; pending wet-lab infiltration trial.",
                validation_level=ValidationLevel.NOT_VALIDATED,
                timestamp=timestamp,
            )
        )

        # Gated release calculation
        has_hard_fail = any(
            r.severity == Severity.HARD_FAIL and r.result == "FAIL"
            for r in evidence
        )
        valid_for_release = not has_hard_fail

        return Variant(
            variant_id=variant_id,
            parent_variant_id=parent_variant_id,
            design_intent=design_intent,
            interventions=interventions,
            sequence=sequence,
            protein_sequence=protein_sequence,
            gc_percent=gc_val,
            cai_score=cai_val,
            policy_version=self.policy_name,
            valid_for_release=valid_for_release,
            validation_records=evidence,
        )

    def _generate_naive_sp_dna(self, sp_aa: str) -> str:
        """Deterministic naive codon assignment for signal peptide if no seed given."""
        naive_codons = {
            "M": "ATG", "A": "GCA", "C": "TGT", "D": "GAT", "E": "GAA",
            "F": "TTT", "G": "GGA", "H": "CAT", "I": "ATT", "K": "AAA",
            "L": "TTA", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "AGA",
            "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
        }
        return "".join(naive_codons.get(aa, "GCA") for aa in sp_aa)

    def _allocate_mature_cds(
        self,
        mature_aa: str,
        watermark: bool,
    ) -> str:
        """
        Deterministic codon allocation for mature chain using structured CodonPolicy.
        """
        allocated_codons = self.policy.allocate_sequence(mature_aa)

        if not watermark:
            return "".join(allocated_codons)

        # Watermark replacement map: systematic GC-elevating synonymous substitutions
        wm_replace = {
            "A": "GCC", "C": "TGC", "D": "GAC", "E": "GAG", "F": "TTC",
            "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAG", "L": "CTC",
            "N": "AAC", "P": "CCC", "Q": "CAG", "R": "AGG", "S": "AGC",
            "T": "ACC", "V": "GTC", "W": "TGG", "Y": "TAC", "M": "ATG",
        }

        watermarked_codons = []
        for aa, orig_codon in zip(mature_aa, allocated_codons):
            watermarked_codons.append(wm_replace.get(aa, orig_codon))

        return "".join(watermarked_codons)

    def _repair_type_iis_sites(self, dna: str, frozen_boundary_nt: int) -> str:
        """
        Deterministic synonymous auto-repair for Type IIS restriction sites (BsaI, BpiI, BsmBI).
        Mutates only codons outside the frozen boundary, preserving 100% amino acid identity.
        """
        dna = dna.upper()
        max_passes = 40

        for _ in range(max_passes):
            # Scan for any match across all patterns
            found_hit = False
            for check_id, (fwd_pat, rev_pat) in TYPE_IIS_PATTERNS.items():
                for pat in (fwd_pat, rev_pat):
                    for match in pat.finditer(dna):
                        start, end = match.start(), match.end()
                        # Codons overlapping this 6nt recognition site
                        c_start = start // 3
                        c_end = (end - 1) // 3
                        # Candidate editable codon positions (must be strictly >= frozen_boundary_nt)
                        editable = [
                            c_idx for c_idx in range(c_start, c_end + 1)
                            if c_idx * 3 >= frozen_boundary_nt and (c_idx * 3 + 3) <= len(dna)
                        ]
                        if not editable:
                            # Hit is within immutable frozen boundary
                            continue

                        # Try to synonymously mutate one of the editable codons
                        repaired = False
                        for c_idx in editable:
                            curr_codon = dna[c_idx * 3 : c_idx * 3 + 3]
                            aa = GENETIC_CODE.get(curr_codon)
                            if not aa or aa == "*":
                                continue

                            # Preferred synonyms from policy, then genetic code
                            allowed = self.policy.allowed_codons.get(aa, [])
                            synonyms = [c for c in allowed if c != curr_codon]
                            if not synonyms:
                                synonyms = [c for c, a in GENETIC_CODE.items() if a == aa and c != curr_codon]
                            synonyms.sort()  # Deterministic order

                            for syn in synonyms:
                                candidate_dna = dna[:c_idx * 3] + syn + dna[c_idx * 3 + 3:]
                                # Verify site was destroyed and no new site created across local window
                                win_start = max(0, (c_idx - 3) * 3)
                                win_end = min(len(candidate_dna), (c_idx + 4) * 3)
                                local_sub = candidate_dna[win_start:win_end]

                                has_any_site = any(
                                    f_p.search(local_sub) or r_p.search(local_sub)
                                    for f_p, r_p in TYPE_IIS_PATTERNS.values()
                                )
                                if not has_any_site:
                                    dna = candidate_dna
                                    found_hit = True
                                    repaired = True
                                    break
                            if repaired:
                                break
                        if repaired:
                            break
                    if found_hit:
                        break
                if found_hit:
                    break
            if not found_hit:
                # All sites cleared or immutable
                break

        return dna

    def _remediate_splice_sites_locally(self, dna: str, frozen_boundary_nt: int) -> str:
        """
        Motif-local synonymous remediation for cryptic donor splice motifs (AGGT[AG]).
        Mutates ONLY the specific codons overlapping flagged motifs in the mature region,
        preserving global codon balance and amino acid identity.
        """
        dna = dna.upper()
        max_passes = 30
        donor_pattern = re.compile(r"AGGT[AG]", re.IGNORECASE)

        for _ in range(max_passes):
            match = donor_pattern.search(dna)
            if not match:
                break

            start, end = match.start(), match.end()
            c_start = start // 3
            c_end = (end - 1) // 3
            editable = [
                c_idx for c_idx in range(c_start, c_end + 1)
                if c_idx * 3 >= frozen_boundary_nt and (c_idx * 3 + 3) <= len(dna)
            ]
            if not editable:
                # Within frozen boundary; cannot mutate
                break

            repaired = False
            for c_idx in editable:
                curr_codon = dna[c_idx * 3 : c_idx * 3 + 3]
                aa = GENETIC_CODE.get(curr_codon)
                if not aa or aa == "*":
                    continue

                allowed = self.policy.allowed_codons.get(aa, [])
                synonyms = [c for c in allowed if c != curr_codon]
                if not synonyms:
                    synonyms = [c for c, a in GENETIC_CODE.items() if a == aa and c != curr_codon]
                synonyms.sort()

                for syn in synonyms:
                    candidate_dna = dna[:c_idx * 3] + syn + dna[c_idx * 3 + 3:]
                    # Check if donor motif is cleared in local window
                    win_start = max(0, (c_idx - 3) * 3)
                    win_end = min(len(candidate_dna), (c_idx + 4) * 3)
                    local_sub = candidate_dna[win_start:win_end]

                    if not donor_pattern.search(local_sub):
                        # Ensure we did not introduce Type IIS site
                        has_iis = any(
                            f_p.search(local_sub) or r_p.search(local_sub)
                            for f_p, r_p in TYPE_IIS_PATTERNS.values()
                        )
                        if not has_iis:
                            dna = candidate_dna
                            repaired = True
                            break
                if repaired:
                    break

            if not repaired:
                break

        return dna
