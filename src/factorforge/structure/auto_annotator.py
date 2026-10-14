"""Reviewable 3D Structure Pause/Linker Recommendation Adapter.

Implements Job 395 P3 requirements:
- Exact full-chain mapping invariant: validates 100% exact residue identity between structure and protein.
- Provenance & confidence discernment: does not treat experimental B-factors as AlphaFold pLDDT.
- Heuristic candidate regions: identifies confidence observations for human review.
- Strict human-in-the-loop contract: outputs reviewed=False, requiring explicit confirmation before use.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from io import StringIO
import math
from typing import Any, Dict, List, Optional
import warnings

from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1


RECOMMENDER_VERSION = "1.0.0"
LIMITATIONS_STATEMENT = (
    "Heuristic structural annotations are recommendations for human review only. "
    "B-factor and low-pLDDT observations do not assign confirmed loops or calibrate translation velocity "
    "or guarantee co-translational folding accuracy. Explicit human review is required."
)


@dataclass
class CandidateRegion:
    start_1based: int
    end_1based: int
    length: int
    region_type: str  # e.g., "flexible_loop", "interdomain_linker", "low_confidence_coil"
    reason: str
    heuristic_basis: str


@dataclass
class StructureRecommendationReport:
    """Full recommendation payload for human inspection."""

    protein_length: int
    chain_id: str
    model_number: int
    format: str
    exact_sequence_match: bool
    is_alphafold_provenance: bool
    confidence_basis: str
    candidate_regions: List[CandidateRegion] = field(default_factory=list)
    recommended_positions: List[int] = field(default_factory=list)
    coverage: float = 0.0
    tool_version: str = RECOMMENDER_VERSION
    limitations: str = LIMITATIONS_STATEMENT
    status: str = "AVAILABLE"
    reason: Optional[str] = None
    provenance_basis: str = "unknown"
    position_mapping: List[Dict[str, Any]] = field(default_factory=list)
    reviewed: bool = False  # Hard invariant: Always False until explicit human confirmation

    def to_dict(self) -> Dict[str, Any]:
        return {
            "protein_length": self.protein_length,
            "chain_id": self.chain_id,
            "model_number": self.model_number,
            "format": self.format,
            "exact_sequence_match": self.exact_sequence_match,
            "is_alphafold_provenance": self.is_alphafold_provenance,
            "confidence_basis": self.confidence_basis,
            "candidate_regions": [asdict(r) for r in self.candidate_regions],
            "recommended_positions": self.recommended_positions,
            "coverage": self.coverage,
            "tool_version": self.tool_version,
            "limitations": self.limitations,
            "reviewed": self.reviewed,
            "status": self.status,
            "reason": self.reason,
            "provenance_basis": self.provenance_basis,
            "position_mapping": self.position_mapping,
        }


def recommend_structure_pause_regions(
    protein: str,
    structure_content: str,
    format: str = "pdb",
    chain_id: str = "A",
    model_number: int = 1,
    is_alphafold: Optional[bool] = None,
    min_region_len: int = 3,
    max_region_len: int = 15,
) -> StructureRecommendationReport:
    """Analyze uploaded structure and recommend candidate pause/linker regions for human review."""
    if not isinstance(protein, str) or not isinstance(format, str):
        raise ValueError("Protein and structure format must be strings")
    clean_prot = "".join(protein.upper().split()).rstrip("*")
    if set(clean_prot) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("Protein must contain unambiguous amino acids")
    if not 3 <= len(clean_prot) <= 300:
        raise ValueError(
            "Structure recommendation supports sequences between 3 and 300 amino acids."
        )

    fmt = format.lower()
    if fmt not in {"pdb", "mmcif"}:
        raise ValueError("Structure format must be 'pdb' or 'mmcif'.")

    if (
        not isinstance(structure_content, str)
        or not structure_content
        or len(structure_content.encode("utf-8")) > 250000
    ):
        raise ValueError("Structure must be a PDB/mmCIF text file within 250 KB")
    if (
        type(model_number) is not int
        or model_number < 1
        or not isinstance(chain_id, str)
        or not 1 <= len(chain_id) <= 8
    ):
        raise ValueError("Invalid explicit chain/model")
    if (
        type(min_region_len) is not int
        or type(max_region_len) is not int
        or not 1 <= min_region_len <= max_region_len <= 300
    ):
        raise ValueError("Invalid region length bounds")
    if is_alphafold is not None and type(is_alphafold) is not bool:
        raise ValueError("AlphaFold provenance flag must be boolean or null")
    # Parse structure
    try:
        parser = (
            PDBParser(PERMISSIVE=False, QUIET=True) if fmt == "pdb" else MMCIFParser(QUIET=True)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure("rec_input", StringIO(structure_content))
        models = list(structure.get_models())
        if model_number > len(models) or model_number < 1:
            raise ValueError(
                f"Model {model_number} not found (structure has {len(models)} models)."
            )
        chain = models[model_number - 1][chain_id]
        residues = [r for r in chain.get_residues() if r.id[0] == " "]
        if any(r.is_disordered() == 2 for r in residues):
            raise ValueError("Ambiguous residue identities")
        observed_seq = "".join(seq1(r.resname, undef_code="X") for r in residues)
    except Exception:
        raise ValueError("Cannot read the selected structure chain/model") from None

    # Invariant: Exact full-chain sequence match
    if observed_seq != clean_prot:
        raise ValueError(
            f"Selected structure chain '{chain_id}' sequence does not match input protein "
            f"(observed length {len(observed_seq)}, input length {len(clean_prot)}). "
            f"Exact full-chain match is required."
        )

    # Discern provenance & B-factors vs pLDDT
    is_af = is_alphafold
    if is_af is None:
        # Check if structure header indicates AlphaFold prediction
        header_text = structure_content[:2000].lower()
        is_af = "alphafold" in header_text
    provenance_basis = (
        "user_declared" if is_alphafold is not None else "file_header_claim" if is_af else "unknown"
    )

    def unavailable(reason):
        return StructureRecommendationReport(
            protein_length=len(clean_prot),
            chain_id=chain_id,
            model_number=model_number,
            format=fmt,
            exact_sequence_match=True,
            is_alphafold_provenance=bool(is_af),
            confidence_basis="not_available",
            status="NOT_AVAILABLE",
            reason=reason,
            provenance_basis=provenance_basis,
        )

    if any("CA" not in r for r in residues):
        return unavailable("CA observations are missing; enter reviewed positions manually")
    b_factors = [float(r["CA"].get_bfactor()) for r in residues]
    if any(not math.isfinite(v) or v < 0 or (is_af and v > 100) for v in b_factors):
        return unavailable("Confidence observations are invalid or outside their declared range")
    if provenance_basis == "unknown":
        return unavailable("Declare the confidence source before using B-factor observations")

    if is_af:
        confidence_basis = (
            "AlphaFold per-residue pLDDT (predicted local distance difference test, 0-100)"
        )
    else:
        confidence_basis = "Experimental crystallographic/NMR B-factors (not calibrated pLDDT)"

    # Identify candidate pause/loop regions
    # Heuristic 1: If AlphaFold, flexible/disordered regions correspond to pLDDT < 70
    # Heuristic 2: Relative B-factor observations; no DSSP or geometric coil assignment
    candidate_regions: List[CandidateRegion] = []

    # Analyze flexibility / loop regions
    flexible_indices = []
    for idx, (res, bf) in enumerate(zip(residues, b_factors), start=1):
        if is_af:
            # Low pLDDT is suggestive of flexible loops or disorder
            if bf < 70.0:
                flexible_indices.append(idx)
        else:
            # For experimental structures, residues with high relative B-factor compared to mean
            mean_b = sum(b_factors) / max(1, len(b_factors))
            if bf > mean_b * 1.2 and bf > 20.0:
                flexible_indices.append(idx)

    # Segment flexible indices into contiguous runs
    runs: List[List[int]] = []
    current_run: List[int] = []
    for pos in flexible_indices:
        if not current_run:
            current_run.append(pos)
        elif pos == current_run[-1] + 1:
            current_run.append(pos)
        else:
            runs.append(current_run)
            current_run = [pos]
    if current_run:
        runs.append(current_run)

    recommended_positions: List[int] = []

    for run in runs:
        run_len = len(run)
        if min_region_len <= run_len <= max_region_len:
            region = CandidateRegion(
                start_1based=run[0],
                end_1based=run[-1],
                length=run_len,
                region_type="low_confidence_region_candidate"
                if is_af
                else "relative_bfactor_region_candidate",
                reason=f"Contiguous observation-based segment of {run_len} residues",
                heuristic_basis=f"{confidence_basis}; run spanning {run[0]}..{run[-1]}",
            )
            candidate_regions.append(region)
            recommended_positions.extend(run)
        elif run_len > max_region_len:
            # Subdivide overly long runs into bounded recommended segments
            for sub_start in range(run[0], run[-1] + 1, max_region_len):
                sub_end = min(sub_start + max_region_len - 1, run[-1])
                sub_len = sub_end - sub_start + 1
                if sub_len >= min_region_len:
                    region = CandidateRegion(
                        start_1based=sub_start,
                        end_1based=sub_end,
                        length=sub_len,
                        region_type="segmented_observation_region_candidate",
                        reason=f"Bounded sub-segment of {sub_len} residues from an extended observation region",
                        heuristic_basis=f"{confidence_basis}; window spanning {sub_start}..{sub_end}",
                    )
                    candidate_regions.append(region)
                    recommended_positions.extend(range(sub_start, sub_end + 1))

    unique_positions = sorted(set(recommended_positions))
    coverage = round(len(unique_positions) / max(1, len(clean_prot)), 4)

    return StructureRecommendationReport(
        protein_length=len(clean_prot),
        chain_id=chain_id,
        model_number=model_number,
        format=fmt,
        exact_sequence_match=True,
        is_alphafold_provenance=bool(is_af),
        confidence_basis=confidence_basis,
        candidate_regions=candidate_regions,
        recommended_positions=unique_positions,
        coverage=coverage,
        tool_version=RECOMMENDER_VERSION,
        limitations=LIMITATIONS_STATEMENT,
        status="AVAILABLE" if candidate_regions else "NO_CANDIDATES",
        reason=None
        if candidate_regions
        else "No supported region detected; enter reviewed positions manually",
        provenance_basis=provenance_basis,
        position_mapping=[
            {
                "protein_position": i,
                "auth_seq_id": r.id[1],
                "insertion_code": r.id[2].strip(),
                "chain_id": chain_id,
                "model_number": model_number,
            }
            for i, r in enumerate(residues, 1)
        ],
        reviewed=False,  # Human review invariant
    )
