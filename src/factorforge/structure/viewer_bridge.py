"""Biomolecular Structure Viewer MCP Bridge for FactorForge.

Translates FactorForge structural recommendation payloads (auto_annotator,
esm_annotator, uniprot_sidecar_fetcher, and interpro_sidecar_fetcher) into
standardized call parameters for the `science_biomolecular_structure_viewer_mcp` sidecar tools.

Key Invariants:
- Distinguishes verified 3D structure residue coordinates (auth_seq_id) from
  unverified 1-based sequence positions.
- Strict rejection: If structure sequence does not correspond to the query protein,
  mapping is rejected and never falls back to blind sequential numbering.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence
import warnings
from io import StringIO

try:
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.SeqUtils import seq1
    _BIOPYTHON_AVAILABLE = True
except ImportError:
    _BIOPYTHON_AVAILABLE = False


COLOR_PALETTE = {
    "predicted_flexible_hinge": "#ff9900",       # Bright orange
    "uniprot_signal_peptide": "#ffd700",         # Gold / yellow
    "uniprot_transmembrane_domain": "#00d4ff",   # Vivid cyan
    "uniprot_disulfide_bond": "#00e676",         # Bright green
    "interpro_domain": "#9c27b0",                # Purple
    "interdomain_boundary": "#e91e63",           # Deep rose
    "low_confidence_region_candidate": "#ff007f",# Hot pink / magenta
    "interdomain_linker": "#ffaa00",
    "flexible_loop": "#ff7700",
    "candidate_region": "#ff00aa",
}


def build_structure_position_mapping(
    protein_sequence: str,
    structure_content: str,
    format: str = "cif",
    chain_id: str = "A",
    model_number: int = 1,
) -> List[Dict[str, Any]]:
    """Build exact 1-based sequence position to structure residue (auth_seq_id) mapping.

    Strict Invariants:
    - Never uses a blind sequential fallback.
    - If the structure sequence does not correspond to the query protein sequence,
      raises ValueError so unverified coordinates are never falsely marked as verified.
    """
    if not _BIOPYTHON_AVAILABLE:
        raise RuntimeError("BioPython is required to parse 3D structure coordinate mapping.")

    clean_prot = "".join(protein_sequence.upper().split()).rstrip("*")
    if not clean_prot:
        raise ValueError("Cannot map coordinates for empty protein sequence.")

    fmt = format.lower().lstrip(".")
    if fmt in {"cif", "mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif fmt == "pdb":
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format '{format}'. Use 'cif' or 'pdb'.")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("struct_map", StringIO(structure_content))

    models = list(structure.get_models())
    if model_number > len(models) or model_number < 1:
        raise ValueError(f"Model {model_number} not found (structure has {len(models)} models).")

    model = models[model_number - 1]
    if chain_id not in model:
        available_chains = [c.id for c in model.get_chains()]
        raise ValueError(f"Chain '{chain_id}' not found in structure. Available chains: {available_chains}")

    chain = model[chain_id]
    residues = [r for r in chain.get_residues() if r.id[0] == " "]
    if not residues:
        raise ValueError(f"Structure chain '{chain_id}' has no valid amino acid residues.")

    observed_seq = "".join(seq1(r.resname, undef_code="X") for r in residues)

    mapping: List[Dict[str, Any]] = []

    # Case 1: Exact full-length sequence match
    if observed_seq == clean_prot:
        for i, r in enumerate(residues, 1):
            mapping.append({
                "protein_position": i,
                "auth_seq_id": r.id[1],
                "insertion_code": r.id[2].strip(),
                "chain_id": chain_id,
                "model_number": model_number,
                "resname": r.resname,
            })
        return mapping

    # Case 2: Exact substring containment (clean_prot in observed_seq)
    # E.g. structure contains expression/purification tags (e.g. His-tag) or fusion partner
    offset = observed_seq.find(clean_prot)
    if offset != -1:
        matched_residues = residues[offset : offset + len(clean_prot)]
        for i, r in enumerate(matched_residues, 1):
            mapping.append({
                "protein_position": i,
                "auth_seq_id": r.id[1],
                "insertion_code": r.id[2].strip(),
                "chain_id": chain_id,
                "model_number": model_number,
                "resname": r.resname,
            })
        return mapping

    # Case 3: Structure is an exact sub-region / domain of the full input protein
    rev_offset = clean_prot.find(observed_seq)
    if rev_offset != -1:
        for j, r in enumerate(residues):
            mapping.append({
                "protein_position": rev_offset + 1 + j,
                "auth_seq_id": r.id[1],
                "insertion_code": r.id[2].strip(),
                "chain_id": chain_id,
                "model_number": model_number,
                "resname": r.resname,
            })
        return mapping

    # Case 4: High-identity pairwise alignment (e.g. point mutations or terminal clipping >= 85% identity)
    try:
        from Bio.Align import PairwiseAligner
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5

        shorter_len = min(len(clean_prot), len(observed_seq))
        score = aligner.score(clean_prot, observed_seq)

        if score >= 0.85 * shorter_len:
            # Retrieve single optimal alignment
            for best_aln in aligner.align(clean_prot, observed_seq):
                aligned_a, aligned_b = best_aln[0], best_aln[1]
                prot_idx = 0
                obs_idx = 0
                for a, b in zip(aligned_a, aligned_b):
                    if a != "-":
                        prot_idx += 1
                    if b != "-":
                        obs_idx += 1
                    if a != "-" and b != "-" and a == b:
                        r = residues[obs_idx - 1]
                        mapping.append({
                            "protein_position": prot_idx,
                            "auth_seq_id": r.id[1],
                            "insertion_code": r.id[2].strip(),
                            "chain_id": chain_id,
                            "model_number": model_number,
                            "resname": r.resname,
                        })
                break
            if mapping:
                return mapping
    except Exception:
        pass

    # Strict rejection: never fallback to blind sequential numbering
    raise ValueError(
        f"Sequence mismatch: structure chain '{chain_id}' sequence (length {len(observed_seq)}) "
        f"does not correspond to query protein sequence (length {len(clean_prot)}). "
        "Verified 3D coordinate mapping rejected due to lack of sequence correspondence."
    )


def prepare_viewer_annotations(
    report_data: Any,
    auth_asym_id: str = "A",
    position_mapping: Optional[Sequence[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Convert candidate regions from recommendation reports into Mol* viewer annotation payloads.

    Validates coordinates against `position_mapping` when available to accurately
    resolve `auth_seq_id` and prevent display errors caused by sequence offsets or truncations.

    Returns a list of dicts suitable for passing directly to `annotate_region` tool.
    """
    candidate_regions = []
    mapping_source = position_mapping

    # Handle object or dict
    if hasattr(report_data, "candidate_regions"):
        candidate_regions = report_data.candidate_regions
        if mapping_source is None and hasattr(report_data, "position_mapping"):
            mapping_source = getattr(report_data, "position_mapping", None)
    elif isinstance(report_data, dict):
        candidate_regions = report_data.get("candidate_regions", [])
        if mapping_source is None:
            mapping_source = report_data.get("position_mapping")

    # Index position mapping if provided
    pos_to_auth: Dict[int, Dict[str, Any]] = {}
    if mapping_source:
        for entry in mapping_source:
            if isinstance(entry, dict) and "protein_position" in entry and "auth_seq_id" in entry:
                pos_to_auth[entry["protein_position"]] = entry

    annotations: List[Dict[str, Any]] = []

    for cr in candidate_regions:
        if hasattr(cr, "start_1based"):
            s = cr.start_1based
            e = cr.end_1based
            rtype = cr.region_type
            reason = cr.reason
        elif isinstance(cr, dict):
            s = cr.get("start_1based")
            e = cr.get("end_1based")
            rtype = cr.get("region_type", "candidate_region")
            reason = cr.get("reason", "")
        else:
            continue

        if s is None or e is None:
            continue

        color = COLOR_PALETTE.get(rtype, COLOR_PALETTE["candidate_region"])
        chain = auth_asym_id
        start_auth = int(s)
        end_auth = int(e)
        coord_mode = "ASSUMED_1_TO_1_UNVERIFIED"
        map_status = "UNVERIFIED"

        if pos_to_auth:
            # Resolve actual structure coordinates
            s_entry = pos_to_auth.get(int(s))
            e_entry = pos_to_auth.get(int(e))

            if s_entry is not None and e_entry is not None:
                start_auth = int(s_entry["auth_seq_id"])
                end_auth = int(e_entry["auth_seq_id"])
                chain = s_entry.get("chain_id", chain)
                coord_mode = "VERIFIED_STRUCTURE_AUTH_SEQ_ID"
                map_status = "EXACT"
            else:
                # Region residues partially present or offset
                present_pos = [p for p in range(int(s), int(e) + 1) if p in pos_to_auth]
                if present_pos:
                    start_auth = int(pos_to_auth[min(present_pos)]["auth_seq_id"])
                    end_auth = int(pos_to_auth[max(present_pos)]["auth_seq_id"])
                    chain = pos_to_auth[min(present_pos)].get("chain_id", chain)
                    coord_mode = "VERIFIED_STRUCTURE_AUTH_SEQ_ID"
                    map_status = f"PARTIAL_{len(present_pos)}_OF_{int(e) - int(s) + 1}_RESIDUES"
                else:
                    # Residues entirely missing from structure (e.g. disordered loop or deleted terminus)
                    coord_mode = "UNMAPPED_RESIDUES_MISSING_FROM_STRUCTURE"
                    map_status = "UNMAPPED"

        if coord_mode == "VERIFIED_STRUCTURE_AUTH_SEQ_ID":
            label = f"{rtype.replace('_', ' ').title()}: {reason} [Seq: {s}..{e} -> Struct: {start_auth}..{end_auth}]"
        elif coord_mode == "UNMAPPED_RESIDUES_MISSING_FROM_STRUCTURE":
            label = f"{rtype.replace('_', ' ').title()}: {reason} [Seq: {s}..{e} (Missing in 3D Structure)]"
        else:
            label = f"{rtype.replace('_', ' ').title()}: {reason} [Seq: {s}..{e} (Structure Unverified)]"

        annotations.append({
            "auth_asym_id": chain,
            "start_auth_seq_id": start_auth,
            "end_auth_seq_id": end_auth,
            "color": color,
            "text": label,
            "coordinate_mapping": coord_mode,
            "mapping_status": map_status,
            "sequence_start_1based": int(s),
            "sequence_end_1based": int(e),
        })

    return annotations
