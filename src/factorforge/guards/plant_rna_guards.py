"""Plant RNA Risk Guards and 1DP vs 3DP Truthful Evaluation Adapter.

Implements Job 394 P0 (truthfulness of measurement/checks) & P1 (unified plant RNA risk scanning):
- Detectors:
  * Cryptic splice donor consensus: [AC]AGGT[AG]AGT (canonical eukaryotic / plant consensus)
  * Premature polyadenylation signals: AATAAA | ATTAAA | AATAAT
  * AU-rich elements (ARE): ATTTA (with tandem/cluster proximity detection)
- Strict evidence classification: SEQUENCE_SCAN on sense transcript
- Truthful labeling: regex match does not guarantee in-vivo processing, and 0 matches does not imply zero biological risk.
- Unified comparison: side-by-side audit of 1DP baseline vs 3DP structure-conditioned candidates.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
import re
from typing import Any, Dict, List, Optional

from factorforge.analysis.metrics import STANDARD_GENETIC_CODE
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs


DETECTOR_VERSION = "1.0.0"
SENSE_SCOPE = "sense_transcript"
EVIDENCE_CLASS = "SEQUENCE_SCAN"

TRUTHFULNESS_NOTICE = (
    "A sequence scan hit indicates a motif pattern match on the sense transcript, "
    "not confirmed biological splicing, polyadenylation, or translation failure. "
    "Conversely, zero scan hits does not guarantee the absence of in-vivo RNA processing risk."
)


@dataclass(frozen=True)
class PlantRNAScanFinding:
    """Individual sequence finding with provenance and exact coordinate bounds."""

    rule_id: str
    name: str
    start_0based: int
    end_exclusive: int
    matched_sequence: str
    scope: str = SENSE_SCOPE
    evidence_class: str = EVIDENCE_CLASS
    detector_version: str = DETECTOR_VERSION
    evaluated_status: str = "evaluated"

    def to_dict(self, include_sequence: bool = True) -> Dict[str, Any]:
        d = asdict(self)
        if not include_sequence:
            d["matched_sequence"] = f"[{len(self.matched_sequence)} nt]"
        return d


@dataclass
class PlantRNAScanResult:
    """Consolidated findings for a single candidate sequence."""

    sequence_length: int
    detector_version: str = DETECTOR_VERSION
    findings: List[PlantRNAScanFinding] = field(default_factory=list)
    donor_count: int = 0
    polya_count: int = 0
    are_count: int = 0
    are_clustered: bool = False
    cluster_details: List[str] = field(default_factory=list)
    summary_status: str = "CLEAN_SCAN"
    truthfulness_notice: str = TRUTHFULNESS_NOTICE

    def to_dict(self, include_sequences: bool = False) -> Dict[str, Any]:
        return {
            "sequence_length": self.sequence_length,
            "detector_version": self.detector_version,
            "evidence_class": EVIDENCE_CLASS,
            "scope": SENSE_SCOPE,
            "summary_status": self.summary_status,
            "donor_count": self.donor_count,
            "polya_count": self.polya_count,
            "are_count": self.are_count,
            "are_clustered": self.are_clustered,
            "cluster_details": self.cluster_details,
            "total_findings": len(self.findings),
            "findings": [f.to_dict(include_sequence=include_sequences) for f in self.findings],
            "truthfulness_notice": self.truthfulness_notice,
        }


def scan_plant_rna_risks(sequence: str) -> PlantRNAScanResult:
    """Scan sense transcript for canonical cryptic donor, polyA, and ARE elements."""
    seq_clean = "".join(sequence.upper().split())
    n = len(seq_clean)
    findings: List[PlantRNAScanFinding] = []

    # 1. Cryptic Splice Donor Consensus: [AC]AGGT[AG]AGT
    donor_pattern = re.compile(r"[AC]AGGT[AG]AGT")
    donor_matches = list(donor_pattern.finditer(seq_clean))
    for m in donor_matches:
        findings.append(
            PlantRNAScanFinding(
                rule_id="rna.cryptic_splice.v1",
                name="Cryptic Splice Donor",
                start_0based=m.start(),
                end_exclusive=m.end(),
                matched_sequence=m.group(),
            )
        )

    # 2. Premature Polyadenylation Signals: AATAAA | ATTAAA | AATAAT
    polya_pattern = re.compile(r"(AATAAA|ATTAAA|AATAAT)")
    polya_matches = list(polya_pattern.finditer(seq_clean))
    for m in polya_matches:
        findings.append(
            PlantRNAScanFinding(
                rule_id="rna.polya_motifs.v1",
                name="Premature PolyA Signal",
                start_0based=m.start(),
                end_exclusive=m.end(),
                matched_sequence=m.group(),
            )
        )

    # 3. AU-rich mRNA Destabilizing Elements (ARE): ATTTA
    are_pattern = re.compile(r"ATTTA")
    are_matches = list(are_pattern.finditer(seq_clean))
    for m in are_matches:
        findings.append(
            PlantRNAScanFinding(
                rule_id="rna.au_rich_elements.v1",
                name="AU-rich Element (ARE)",
                start_0based=m.start(),
                end_exclusive=m.end(),
                matched_sequence=m.group(),
            )
        )

    # Check for ARE clustering: >= 2 ARE motifs within a sliding window of 30 nt
    are_clustered = False
    cluster_details: List[str] = []
    are_positions = [m.start() for m in are_matches]
    for i in range(len(are_positions) - 1):
        gap = are_positions[i + 1] - are_positions[i]
        if gap <= 30:
            are_clustered = True
            cluster_details.append(
                f"ARE cluster at positions {are_positions[i]} and {are_positions[i + 1]} (gap: {gap} nt)"
            )

    total_hits = len(findings)
    status = "CLEAN_SCAN" if total_hits == 0 else "ADVISORY_FINDINGS"

    return PlantRNAScanResult(
        sequence_length=n,
        detector_version=DETECTOR_VERSION,
        findings=findings,
        donor_count=len(donor_matches),
        polya_count=len(polya_matches),
        are_count=len(are_matches),
        are_clustered=are_clustered,
        cluster_details=cluster_details,
        summary_status=status,
        truthfulness_notice=TRUTHFULNESS_NOTICE,
    )


def validate_candidate_hard_invariants(
    cds: str,
    protein: str,
    gc_min: float,
    gc_max: float,
    forbidden_motifs: Optional[List[str]] = None,
    left_flank: str = "",
    right_flank: str = "",
) -> Dict[str, Any]:
    """Independent zero-defect hard invariant validator for a generated candidate.

    Verifies:
    - 100% translation identity (no amino acid substitutions, insertions, deletions)
    - Valid open reading frame (length multiple of 3, starts with ATG if protein starts with M)
    - Zero internal stop codons
    - Exact integer GC count inside [gc_min, gc_max]
    - Zero Type IIS restriction sites (forward and reverse complement) across full construct.
    """
    clean_cds = "".join(cds.upper().split())
    clean_prot = "".join(protein.upper().split()).rstrip("*")

    errors: List[str] = []

    # 1. Length & Reading frame
    if len(clean_cds) != len(clean_prot) * 3:
        errors.append(
            f"Length mismatch: CDS has {len(clean_cds)} nt, expected {len(clean_prot) * 3} nt"
        )
    if len(clean_cds) % 3 != 0:
        errors.append(f"Reading frame error: length {len(clean_cds)} is not a multiple of 3")

    # 2. Translation & Internal Stops
    translated = []
    internal_stops = 0
    for i in range(0, len(clean_cds), 3):
        codon = clean_cds[i : i + 3]
        aa = STANDARD_GENETIC_CODE.get(codon, "?")
        translated.append(aa)
        if aa == "*" and i < len(clean_cds) - 3:
            internal_stops += 1
    translated_str = "".join(translated)

    if translated_str != clean_prot:
        errors.append(
            f"Translation mismatch: expected {clean_prot[:10]}... but got {translated_str[:10]}..."
        )
    if internal_stops > 0:
        errors.append(f"Found {internal_stops} internal stop codon(s)")

    # 3. Exact Integer GC Bounds
    gc_count = clean_cds.count("G") + clean_cds.count("C")
    lower_bound = gc_min if gc_min <= 1.0 else gc_min / 100.0
    upper_bound = gc_max if gc_max <= 1.0 else gc_max / 100.0
    min_gc_count = math.ceil(lower_bound * len(clean_cds))
    max_gc_count = math.floor(upper_bound * len(clean_cds))

    gc_feasible = min_gc_count <= gc_count <= max_gc_count
    actual_gc_fraction = round(gc_count / max(1, len(clean_cds)), 4)
    if not gc_feasible:
        errors.append(
            f"GC count {gc_count} ({actual_gc_fraction * 100:.1f}%) is outside "
            f"allowed range [{min_gc_count}, {max_gc_count}] ({lower_bound * 100:.1f}%-{upper_bound * 100:.1f}%)"
        )

    # 4. Forbidden Type IIS Motifs (Sense & Reverse Complement)
    motifs = get_canonical_forbidden_motifs() if forbidden_motifs is None else forbidden_motifs
    construct = (left_flank + clean_cds + right_flank).upper()
    complement_map = str.maketrans("ACGT", "TGCA")
    found_motifs: List[str] = []

    for m in motifs:
        m_upper = m.upper()
        m_rc = m_upper.translate(complement_map)[::-1]
        if m_upper in construct or m_rc in construct:
            found_motifs.append(m_upper)

    if found_motifs:
        errors.append(f"Forbidden restriction site(s) found in construct: {found_motifs}")

    passed = len(errors) == 0

    return {
        "passed": passed,
        "translation_identity_100pct": (translated_str == clean_prot),
        "actual_gc_fraction": actual_gc_fraction,
        "actual_gc_count": gc_count,
        "allowed_gc_range": [min_gc_count, max_gc_count],
        "gc_feasible": gc_feasible,
        "internal_stops": internal_stops,
        "forbidden_motifs_found": found_motifs,
        "errors": errors,
    }


def compare_candidate_plant_rna_risks(
    cds_1dp: str,
    cds_3dp: str,
    protein_name: Optional[str] = None,
) -> Dict[str, Any]:
    """Execute unified plant RNA risk scanning across both 1DP and 3DP candidates."""
    res_1dp = scan_plant_rna_risks(cds_1dp)
    res_3dp = scan_plant_rna_risks(cds_3dp)

    delta_donors = res_3dp.donor_count - res_1dp.donor_count
    delta_polya = res_3dp.polya_count - res_1dp.polya_count
    delta_are = res_3dp.are_count - res_1dp.are_count

    # Comparison summary text table
    sign_donor = "+" if delta_donors > 0 else ""
    sign_polya = "+" if delta_polya > 0 else ""
    sign_are = "+" if delta_are > 0 else ""

    summary_lines = [
        f"Plant RNA Risk Comparison ({protein_name or 'Candidate'}):",
        "  Feature                1DP     3DP     Delta",
        f"  Cryptic Donors:        {res_1dp.donor_count:3d}     {res_3dp.donor_count:3d}     {sign_donor}{delta_donors:d}",
        f"  Premature PolyA:       {res_1dp.polya_count:3d}     {res_3dp.polya_count:3d}     {sign_polya}{delta_polya:d}",
        f"  AU-rich Elements:      {res_1dp.are_count:3d}     {res_3dp.are_count:3d}     {sign_are}{delta_are:d}",
        f"  ARE Clustering:        {str(res_1dp.are_clustered):5s}   {str(res_3dp.are_clustered):5s}   --",
        f"  Overall Scan Status:   {res_1dp.summary_status} vs {res_3dp.summary_status}",
        "",
        f"Note: {TRUTHFULNESS_NOTICE}",
    ]
    summary_table = "\n".join(summary_lines)

    return {
        "protein_name": protein_name,
        "detector_version": DETECTOR_VERSION,
        "1dp_scan": res_1dp.to_dict(include_sequences=False),
        "3dp_scan": res_3dp.to_dict(include_sequences=False),
        "delta_donors": delta_donors,
        "delta_polya": delta_polya,
        "delta_are": delta_are,
        "are_clustered_1dp": res_1dp.are_clustered,
        "are_clustered_3dp": res_3dp.are_clustered,
        "summary_table": summary_table,
        "truthfulness_notice": TRUTHFULNESS_NOTICE,
    }
