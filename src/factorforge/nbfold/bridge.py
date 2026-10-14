"""StructureBridge: Bounded AF2/AF3/PDB/mmCIF intake, CA-confidence observation, and sidecar validation.

Implements PRD-JOB-398-NB-FOLD-v1.1.0 Section 4:
- Size checks before full file reads (FASTA <= 16 KB, structure <= 250 KB, sidecar <= 2 MB, total <= 5 MB).
- Single-record standard 20-AA protein validation (3-300 aa).
- Format/extension normalization (.cif/.mmcif -> mmcif, .pdb -> pdb) and content cross-checking.
- Explicit chain/model selection and 100% exact full-chain sequence identity/length match.
- Strict rejection of modified, hetero, or disordered residues and missing CA atoms.
- Explicit distinction between AF2 residue pLDDT, AF3 CA-atom confidence, and experimental B-factors.
- Optional confidence/PAE sidecar validation with explicit chain/token/atom mapping.
- Reuses factorforge.structure.auto_annotator heuristic without inventing midpoint fallbacks.
"""

from __future__ import annotations

import json
import math
import re
import warnings
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFParser, PDBParser

from factorforge.nbfold.annotations import (
    analyze_backbone_secondary_structure,
    analyze_domain_annotations,
    analyze_uniprot_tm_signal_annotations,
    build_structural_review_summary,
)
from factorforge.nbfold.contracts import (
    ANNOTATION_SCHEMA_VERSION,
    MAX_FASTA_BYTES,
    MAX_MANIFEST_JSON_BYTES,
    MAX_PROTEIN_AA,
    MAX_SIDECAR_BYTES,
    MAX_STRUCTURE_BYTES,
    MAX_TOTAL_INPUT_BYTES,
    MIN_PROTEIN_AA,
    NBFOLD_PIPELINE_VERSION,
    NBFoldError,
    check_file_size_limit,
    compute_annotation_content_sha256,
    sha256_bytes,
    sha256_text,
)
from factorforge.structure.auto_annotator import (
    LIMITATIONS_STATEMENT,
    RECOMMENDER_VERSION,
    recommend_structure_pause_regions,
)

VALID_STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

CANONICAL_3TO1: dict[str, str] = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

VALID_PROVIDER_FAMILIES = frozenset({"af2", "af3", "experimental", "unknown"})


@dataclass(frozen=True)
class ParsedProteinInput:
    header: str
    sequence: str
    protein_sha256: str
    file_sha256: str
    byte_size: int


def parse_protein_fasta_file(fasta_path: Path, *, stage: str = "annotate") -> ParsedProteinInput:
    """Validate size and parse a single-record standard-20-AA protein FASTA file."""
    p = Path(fasta_path)
    byte_size = check_file_size_limit(
        p,
        max_bytes=MAX_FASTA_BYTES,
        label="Protein FASTA",
        empty_error_code="INVALID_FASTA",
        stage=stage,
    )
    raw_bytes = p.read_bytes()
    file_sha256 = sha256_bytes(raw_bytes)
    try:
        raw_text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise NBFoldError(
            "INVALID_FASTA",
            "Protein FASTA file must be valid UTF-8 text.",
            stage=stage,
        ) from exc

    lines = raw_text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    nonempty = [line.strip() for line in lines if line.strip()]
    if not nonempty:
        raise NBFoldError(
            "INVALID_FASTA",
            "Protein FASTA file contains no non-empty lines.",
            stage=stage,
        )

    header_indices = [idx for idx, line in enumerate(nonempty) if line.startswith(">")]
    if len(header_indices) > 1:
        raise NBFoldError(
            "MULTI_RECORD_FASTA",
            "Multiple FASTA records detected. NB-Fold requires a single protein record.",
            stage=stage,
            details={"record_count": len(header_indices)},
        )
    if not header_indices or header_indices[0] != 0:
        raise NBFoldError(
            "INVALID_FASTA",
            "Protein FASTA must start with a '>' header line.",
            stage=stage,
        )

    header = nonempty[0][1:].strip() or "protein_record"
    seq_lines = nonempty[1:]
    if not seq_lines:
        raise NBFoldError(
            "INVALID_SEQUENCE",
            "Protein FASTA record contains no sequence lines.",
            stage=stage,
        )

    compact = re.sub(r"\s+", "", "".join(seq_lines)).upper()
    if not compact:
        raise NBFoldError(
            "INVALID_SEQUENCE",
            "Protein sequence is empty after whitespace normalization.",
            stage=stage,
        )

    if "*" in compact[:-1]:
        raise NBFoldError(
            "INVALID_SEQUENCE",
            "Protein sequence contains an internal stop codon ('*').",
            stage=stage,
        )

    protein = compact.rstrip("*")
    if not protein:
        raise NBFoldError(
            "INVALID_SEQUENCE",
            "Protein sequence must contain standard amino acid residues.",
            stage=stage,
        )

    invalid_chars = sorted(set(protein) - VALID_STANDARD_AA)
    if invalid_chars:
        raise NBFoldError(
            "INVALID_SEQUENCE",
            f"Protein sequence contains non-standard or ambiguous character(s): {', '.join(invalid_chars)}.",
            stage=stage,
            details={"invalid_characters": invalid_chars},
        )

    if not (MIN_PROTEIN_AA <= len(protein) <= MAX_PROTEIN_AA):
        raise NBFoldError(
            "UNSUPPORTED_INPUT_SIZE",
            f"Protein length ({len(protein)} aa) is outside supported Phase A/B bounds "
            f"({MIN_PROTEIN_AA}-{MAX_PROTEIN_AA} aa).",
            stage=stage,
            details={
                "protein_length_aa": len(protein),
                "min_aa": MIN_PROTEIN_AA,
                "max_aa": MAX_PROTEIN_AA,
            },
        )

    return ParsedProteinInput(
        header=header,
        sequence=protein,
        protein_sha256=sha256_text(protein),
        file_sha256=file_sha256,
        byte_size=byte_size,
    )


def resolve_and_validate_structure_format(
    structure_path: Path,
    declared_format: str | None = None,
    *,
    stage: str = "annotate",
) -> tuple[str, str, str, int]:
    """Validate structure file size, extension, declared format, and content consistency.

    Returns (normalized_format, structure_text, structure_sha256, byte_size).
    """
    p = Path(structure_path)
    byte_size = check_file_size_limit(
        p,
        max_bytes=MAX_STRUCTURE_BYTES,
        label="Structure",
        empty_error_code="INVALID_STRUCTURE_FORMAT",
        stage=stage,
    )
    ext = p.suffix.lower()
    if ext == ".pdb":
        ext_format = "pdb"
    elif ext in {".cif", ".mmcif"}:
        ext_format = "mmcif"
    else:
        raise NBFoldError(
            "INVALID_STRUCTURE_FORMAT",
            f"Unsupported structure file extension '{ext}'. Expected .pdb, .cif, or .mmcif.",
            stage=stage,
        )

    if declared_format is not None:
        if not isinstance(declared_format, str) or not declared_format.strip():
            raise NBFoldError(
                "INVALID_STRUCTURE_FORMAT",
                "Structure format must be a non-empty string ('pdb', 'cif', or 'mmcif').",
                stage=stage,
            )
        norm_decl = declared_format.strip().lower()
        if norm_decl in {"cif", "mmcif"}:
            norm_decl = "mmcif"
        elif norm_decl == "pdb":
            norm_decl = "pdb"
        else:
            raise NBFoldError(
                "INVALID_STRUCTURE_FORMAT",
                f"Unsupported declared format '{declared_format}'. Use 'pdb', 'cif', or 'mmcif'.",
                stage=stage,
            )
        if norm_decl != ext_format:
            raise NBFoldError(
                "FORMAT_EXTENSION_MISMATCH",
                f"Declared format '{declared_format}' (normalized '{norm_decl}') conflicts "
                f"with file extension '{ext}' ('{ext_format}').",
                stage=stage,
            )
        fmt = norm_decl
    else:
        fmt = ext_format

    raw_bytes = p.read_bytes()
    structure_sha256 = sha256_bytes(raw_bytes)
    try:
        text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise NBFoldError(
            "INVALID_STRUCTURE_FORMAT",
            "Structure file must be valid UTF-8 text.",
            stage=stage,
        ) from exc

    stripped = text.lstrip()
    has_mmcif_markers = stripped.startswith("data_") or "_atom_site." in text
    has_pdb_atom_lines = any(
        line.startswith(("ATOM  ", "HETATM")) for line in text.splitlines()
    )

    if fmt == "pdb":
        if has_mmcif_markers:
            raise NBFoldError(
                "FORMAT_EXTENSION_MISMATCH",
                "Structure file is declared/named as PDB but contains mmCIF data blocks.",
                stage=stage,
            )
        if not has_pdb_atom_lines:
            raise NBFoldError(
                "INVALID_STRUCTURE_FORMAT",
                "PDB structure file contains no ATOM records.",
                stage=stage,
            )
    else:
        if not has_mmcif_markers:
            if has_pdb_atom_lines:
                raise NBFoldError(
                    "FORMAT_EXTENSION_MISMATCH",
                    "Structure file is declared/named as mmCIF but contains fixed-column PDB records.",
                    stage=stage,
                )
            raise NBFoldError(
                "INVALID_STRUCTURE_FORMAT",
                "mmCIF structure file is missing required data_ / _atom_site records.",
                stage=stage,
            )

    return fmt, text, structure_sha256, byte_size


def _prescan_nonfinite_bfactors(structure_text: str, fmt: str, *, stage: str) -> None:
    """Detect literal NaN/Inf tokens in B-factor columns before Biopython parser fails generically."""
    nonfinite_tokens = {"nan", "+nan", "-nan", "inf", "+inf", "-inf", "infinity", "-infinity"}
    for line in structure_text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if fmt == "pdb" and line.startswith(("ATOM  ", "HETATM")):
            bf_field = line[60:66].strip().lower()
            parts = s.split()
            if bf_field in nonfinite_tokens or any(tok.lower() in nonfinite_tokens for tok in parts):
                raise NBFoldError(
                    "INVALID_CONFIDENCE_VALUE",
                    "Non-finite B-factor/confidence token (NaN/Inf) detected in PDB ATOM record.",
                    stage=stage,
                )
        elif fmt == "mmcif" and s.startswith(("ATOM ", "HETATM ")):
            parts = s.split()
            if any(tok.lower() in nonfinite_tokens for tok in parts):
                raise NBFoldError(
                    "INVALID_CONFIDENCE_VALUE",
                    "Non-finite B-factor/confidence token (NaN/Inf) detected in mmCIF atom record.",
                    stage=stage,
                )


def _extract_mmcif_label_seq_ids(structure_text: str, chain_id: str) -> dict[int, int | None]:
    """Best-effort extraction of (auth_seq_id -> label_seq_id) for CA atoms in mmCIF."""
    mapping: dict[int, int | None] = {}
    lines = structure_text.splitlines()
    headers: list[str] = []
    in_atom_loop = False
    for line in lines:
        s = line.strip()
        if s == "loop_":
            headers = []
            in_atom_loop = False
            continue
        if s.startswith("_atom_site."):
            headers.append(s)
            in_atom_loop = True
            continue
        if in_atom_loop and s and not s.startswith(("_", "#")):
            if s.startswith(("loop_", "data_")):
                in_atom_loop = False
                continue
            parts = s.split()
            if len(parts) < len(headers):
                continue
            row = dict(zip(headers, parts))
            atom_name = row.get("_atom_site.label_atom_id") or row.get("_atom_site.auth_atom_id")
            row_chain = row.get("_atom_site.auth_asym_id") or row.get("_atom_site.label_asym_id")
            if atom_name == "CA" and row_chain == chain_id:
                try:
                    auth_seq = int(
                        row.get("_atom_site.auth_seq_id") or row.get("_atom_site.label_seq_id") or 0
                    )
                    lbl_raw = row.get("_atom_site.label_seq_id")
                    lbl_seq = int(lbl_raw) if lbl_raw and lbl_raw not in {".", "?"} else None
                    mapping[auth_seq] = lbl_seq
                except ValueError:
                    pass
    return mapping


def _resolve_provider_and_source_basis(
    *,
    declared_provider: str | None,
    structure_text: str,
    export_metadata: dict[str, Any] | None,
    stage: str,
) -> tuple[str, str, str | None]:
    """Determine (provider_family, source_basis, export_version) without claiming producer auth."""
    export_version: str | None = None
    meta_provider: str | None = None
    if export_metadata is not None:
        if not isinstance(export_metadata, dict):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "Export metadata sidecar must be a JSON object.",
                stage=stage,
            )
        raw_mp = export_metadata.get("provider_family")
        if raw_mp is not None:
            if not isinstance(raw_mp, str) or raw_mp.strip().lower() not in VALID_PROVIDER_FAMILIES:
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"Invalid provider_family '{raw_mp}' in export metadata.",
                    stage=stage,
                )
            meta_provider = raw_mp.strip().lower()
        raw_ver = export_metadata.get("model_version") or export_metadata.get("export_version")
        if isinstance(raw_ver, str) and raw_ver.strip():
            export_version = raw_ver.strip()

    norm_decl: str | None = None
    if declared_provider is not None:
        if not isinstance(declared_provider, str):
            raise NBFoldError(
                "UNKNOWN_CONFIDENCE_SOURCE",
                "provider_family must be a string ('af2', 'af3', 'experimental', or 'unknown').",
                stage=stage,
            )
        norm_decl = declared_provider.strip().lower()
        if norm_decl not in VALID_PROVIDER_FAMILIES:
            raise NBFoldError(
                "UNKNOWN_CONFIDENCE_SOURCE",
                f"Unsupported provider_family '{declared_provider}'. "
                "Allowed values: af2, af3, experimental, unknown.",
                stage=stage,
            )

    if norm_decl and norm_decl != "unknown":
        if meta_provider and meta_provider != "unknown" and meta_provider != norm_decl:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                f"Declared provider_family '{norm_decl}' conflicts with export metadata '{meta_provider}'.",
                stage=stage,
            )
        basis = "export_metadata" if meta_provider == norm_decl else "user_declared"
        return norm_decl, basis, export_version

    if meta_provider and meta_provider != "unknown":
        return meta_provider, "export_metadata", export_version

    # If user explicitly declared 'unknown', do not guess from header
    if norm_decl == "unknown":
        raise NBFoldError(
            "UNKNOWN_CONFIDENCE_SOURCE",
            "Confidence source is explicitly 'unknown'; declare --provider-family (af2, af3, or experimental) "
            "before using B-factor/confidence observations.",
            stage=stage,
        )

    header_sample = structure_text[:4000].lower()
    if "alphafold 3" in header_sample or "alphafold3" in header_sample or "af3" in header_sample:
        return "af3", "file_header_claim", export_version
    if "alphafold" in header_sample:
        return "af2", "file_header_claim", export_version

    raise NBFoldError(
        "UNKNOWN_CONFIDENCE_SOURCE",
        "Confidence source cannot be determined; declare --provider-family (af2, af3, or experimental).",
        stage=stage,
    )


def _validate_and_map_sidecars(
    *,
    confidence_sidecar_data: Any | None,
    pae_sidecar_data: Any | None,
    chain_id: str,
    protein_length: int,
    ca_bfactors: list[float],
    stage: str,
) -> tuple[dict[str, Any], dict[str, Any], list[int | None]]:
    """Validate optional confidence and PAE sidecars with explicit chain/token/atom mapping."""
    token_indices: list[int | None] = [None] * protein_length
    sidecar_obs: dict[str, Any] = {
        "status": "NOT_AVAILABLE",
        "confidence_sidecar_verified": False,
        "pae_sidecar_verified": False,
        "reason": "No sidecar files provided.",
    }
    domain_rel: dict[str, Any] = {
        "status": "NOT_AVAILABLE",
        "evidence_source": None,
        "used_for_candidate_recommendation": False,
        "domain_segmentation_status": "NOT_AVAILABLE",
        "reason": "No verified PAE sidecar provided.",
    }

    pae_payload = pae_sidecar_data
    if confidence_sidecar_data is not None:
        if not isinstance(confidence_sidecar_data, dict):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "Confidence sidecar must be a JSON object with explicit chain/token/atom mapping.",
                stage=stage,
            )

        # Verify either atom-level mapping or token/residue-level mapping for selected chain
        mapped_ca_values: list[float] | None = None
        mapped_tokens: list[int] | None = None

        if "atom_chain_ids" in confidence_sidecar_data and "atom_plddts" in confidence_sidecar_data:
            atom_chains = confidence_sidecar_data.get("atom_chain_ids")
            atom_names = confidence_sidecar_data.get("atom_names")
            atom_plddts = confidence_sidecar_data.get("atom_plddts")
            if (
                not isinstance(atom_chains, list)
                or not isinstance(atom_names, list)
                or not isinstance(atom_plddts, list)
                or not (len(atom_chains) == len(atom_names) == len(atom_plddts))
            ):
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    "Atom-level confidence sidecar arrays (atom_chain_ids, atom_names, atom_plddts) "
                    "must be equal-length lists.",
                    stage=stage,
                )
            ca_vals: list[float] = []
            tok_idxs: list[int] = []
            for idx, (ch, aname, val) in enumerate(zip(atom_chains, atom_names, atom_plddts)):
                if str(ch) == chain_id and str(aname).strip() == "CA":
                    if (
                        isinstance(val, bool)
                        or not isinstance(val, (int, float))
                        or not math.isfinite(float(val))
                        or not (0.0 <= float(val) <= 100.0)
                    ):
                        raise NBFoldError(
                            "SIDECAR_MISMATCH",
                            "Atom-level sidecar contains non-finite or out-of-range [0, 100] CA confidence.",
                            stage=stage,
                        )
                    ca_vals.append(float(val))
                    tok_idxs.append(idx)
            mapped_ca_values = ca_vals
            mapped_tokens = tok_idxs

        elif "token_chain_ids" in confidence_sidecar_data and (
            "ca_atom_plddts" in confidence_sidecar_data
            or "residue_plddts" in confidence_sidecar_data
            or "plddt" in confidence_sidecar_data
        ):
            tok_chains = confidence_sidecar_data.get("token_chain_ids")
            tok_res_ids = confidence_sidecar_data.get("token_res_ids")
            raw_vals = (
                confidence_sidecar_data.get("ca_atom_plddts")
                if "ca_atom_plddts" in confidence_sidecar_data
                else confidence_sidecar_data.get("residue_plddts")
                if "residue_plddts" in confidence_sidecar_data
                else confidence_sidecar_data.get("plddt")
            )
            if (
                not isinstance(tok_chains, list)
                or not isinstance(tok_res_ids, list)
                or not isinstance(raw_vals, list)
                or not (len(tok_chains) == len(tok_res_ids) == len(raw_vals))
            ):
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    "Token-level confidence sidecar requires equal-length token_chain_ids, "
                    "token_res_ids, and confidence value arrays.",
                    stage=stage,
                )
            ca_vals = []
            tok_idxs = []
            chain_res_ids: list[int] = []
            for idx, (ch, res_id, val) in enumerate(zip(tok_chains, tok_res_ids, raw_vals)):
                if str(ch) == chain_id:
                    if (
                        isinstance(val, bool)
                        or not isinstance(val, (int, float))
                        or not math.isfinite(float(val))
                        or not (0.0 <= float(val) <= 100.0)
                        or type(res_id) is not int
                    ):
                        raise NBFoldError(
                            "SIDECAR_MISMATCH",
                            "Token-level sidecar contains invalid residue ID or non-finite/out-of-range confidence.",
                            stage=stage,
                        )
                    ca_vals.append(float(val))
                    tok_idxs.append(idx)
                    chain_res_ids.append(res_id)
            if chain_res_ids != list(range(1, protein_length + 1)):
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"Sidecar token_res_ids for chain '{chain_id}' do not form an exact 1..{protein_length} mapping.",
                    stage=stage,
                )
            mapped_ca_values = ca_vals
            mapped_tokens = tok_idxs
        else:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "Confidence sidecar lacks explicit chain/token/atom mapping fields "
                "(expected token_chain_ids+token_res_ids or atom_chain_ids+atom_names).",
                stage=stage,
            )

        if len(mapped_ca_values) != protein_length:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                f"Confidence sidecar mapped {len(mapped_ca_values)} residues for chain '{chain_id}', "
                f"expected {protein_length}.",
                stage=stage,
            )

        for pos_idx, (struct_bf, sidecar_bf) in enumerate(
            zip(ca_bfactors, mapped_ca_values), start=1
        ):
            if abs(struct_bf - sidecar_bf) > 0.5:
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"Sidecar confidence at position {pos_idx} ({sidecar_bf:.2f}) conflicts "
                    f"with structure CA B-factor ({struct_bf:.2f}).",
                    stage=stage,
                )

        token_indices = [int(i) for i in mapped_tokens]
        sidecar_obs = {
            "status": "AVAILABLE",
            "confidence_sidecar_verified": True,
            "pae_sidecar_verified": False,
            "mapped_chain_id": chain_id,
            "mapped_residue_count": protein_length,
            "reason": None,
        }
        if pae_payload is None and "pae" in confidence_sidecar_data:
            pae_payload = confidence_sidecar_data

    if pae_payload is not None:
        if not isinstance(pae_payload, dict) or "pae" not in pae_payload:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "PAE sidecar must be a JSON object containing a 'pae' matrix and token_chain_ids.",
                stage=stage,
            )
        pae_matrix = pae_payload.get("pae")
        tok_chains = pae_payload.get("token_chain_ids")
        tok_res_ids = pae_payload.get("token_res_ids")
        if (
            not isinstance(pae_matrix, list)
            or not isinstance(tok_chains, list)
            or not isinstance(tok_res_ids, list)
            or len(tok_chains) != len(tok_res_ids)
            or len(pae_matrix) != len(tok_chains)
        ):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "PAE sidecar matrix dimensions must match token_chain_ids and token_res_ids length.",
                stage=stage,
            )
        n_total = len(tok_chains)
        for row in pae_matrix:
            if not isinstance(row, list) or len(row) != n_total:
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    "PAE matrix must be square (N_tokens x N_tokens).",
                    stage=stage,
                )

        chain_token_indices = [idx for idx, ch in enumerate(tok_chains) if str(ch) == chain_id]
        chain_res = [tok_res_ids[idx] for idx in chain_token_indices]
        if len(chain_token_indices) != protein_length or chain_res != list(
            range(1, protein_length + 1)
        ):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                f"PAE sidecar tokens for chain '{chain_id}' do not match protein length {protein_length}.",
                stage=stage,
            )

        sub_values: list[float] = []
        chain_pae_matrix: list[list[float]] = []
        for i_idx in chain_token_indices:
            row_vals: list[float] = []
            for j_idx in chain_token_indices:
                val = pae_matrix[i_idx][j_idx]
                if (
                    isinstance(val, bool)
                    or not isinstance(val, (int, float))
                    or not math.isfinite(float(val))
                    or float(val) < 0.0
                ):
                    raise NBFoldError(
                        "SIDECAR_MISMATCH",
                        "PAE matrix contains non-finite or negative values.",
                        stage=stage,
                    )
                fval = float(val)
                sub_values.append(fval)
                row_vals.append(fval)
            chain_pae_matrix.append(row_vals)

        if all(t is None for t in token_indices):
            token_indices = [int(i) for i in chain_token_indices]

        mean_pae = round(sum(sub_values) / max(1, len(sub_values)), 3)
        max_pae = round(max(sub_values) if sub_values else 0.0, 3)
        min_pae = round(min(sub_values) if sub_values else 0.0, 3)

        sidecar_obs["status"] = "AVAILABLE"
        sidecar_obs["pae_sidecar_verified"] = True
        sidecar_obs["mapped_chain_id"] = chain_id
        sidecar_obs["mapped_residue_count"] = protein_length
        sidecar_obs["reason"] = None

        domain_rel = {
            "status": "AVAILABLE",
            "evidence_source": "verified_pae_sidecar",
            "chain_id": chain_id,
            "chain_token_count": protein_length,
            "matrix_shape": [protein_length, protein_length],
            "mean_pae_angstrom": mean_pae,
            "min_pae_angstrom": min_pae,
            "max_pae_angstrom": max_pae,
            "used_for_candidate_recommendation": False,
            "domain_segmentation_status": "NOT_AVAILABLE",
            "note": (
                "PAE matrix summary is recorded for human inspection only. "
                "PAE domain segmentation is never used to overwrite auto_annotator candidate positions."
            ),
        }
        return sidecar_obs, domain_rel, token_indices, chain_pae_matrix

    return sidecar_obs, domain_rel, token_indices, None


def build_structure_annotation(
    *,
    protein_fasta_path: Path,
    structure_file_path: Path,
    declared_format: str | None = None,
    chain_id: str = "A",
    model_number: int = 1,
    provider_family: str | None = None,
    confidence_sidecar_path: Path | None = None,
    pae_sidecar_path: Path | None = None,
    export_metadata_path: Path | None = None,
    interpro_sidecar_path: Path | None = None,
    uniprot_sidecar_path: Path | None = None,
    min_region_len: int = 3,
    max_region_len: int = 15,
    default_config_digest: str | None = None,
    resolved_config_preview: dict[str, Any] | None = None,
    stage: str = "annotate",
) -> dict[str, Any]:
    """Execute Phase A/C StructureBridge intake, exact residue mapping, and structural annotation."""
    if not isinstance(chain_id, str) or not (1 <= len(chain_id.strip()) <= 8):
        raise NBFoldError(
            "CHAIN_NOT_FOUND",
            "An explicit structure chain ID (1-8 characters) is required.",
            stage=stage,
        )
    chain_id = chain_id.strip()
    if type(model_number) is not int or model_number < 1:
        raise NBFoldError(
            "MODEL_NOT_FOUND",
            "Model number must be a positive integer (1-based file model order).",
            stage=stage,
        )

    # 1. Pre-check all file sizes BEFORE reading full contents
    total_bytes = 0
    fasta_size = check_file_size_limit(
        Path(protein_fasta_path),
        max_bytes=MAX_FASTA_BYTES,
        label="Protein FASTA",
        empty_error_code="INVALID_FASTA",
        stage=stage,
    )
    total_bytes += fasta_size

    struct_size = check_file_size_limit(
        Path(structure_file_path),
        max_bytes=MAX_STRUCTURE_BYTES,
        label="Structure",
        empty_error_code="INVALID_STRUCTURE_FORMAT",
        stage=stage,
    )
    total_bytes += struct_size

    sidecar_sha256s: dict[str, str | None] = {
        "confidence_sidecar": None,
        "pae_sidecar": None,
        "export_metadata": None,
        "interpro_sidecar": None,
        "uniprot_sidecar": None,
    }
    confidence_sidecar_data: Any | None = None
    pae_sidecar_data: Any | None = None
    export_metadata_data: Any | None = None
    interpro_sidecar_data: Any | None = None
    uniprot_sidecar_data: Any | None = None

    for key, spath, max_b in (
        ("confidence_sidecar", confidence_sidecar_path, MAX_SIDECAR_BYTES),
        ("pae_sidecar", pae_sidecar_path, MAX_SIDECAR_BYTES),
        ("export_metadata", export_metadata_path, MAX_MANIFEST_JSON_BYTES),
        ("interpro_sidecar", interpro_sidecar_path, MAX_SIDECAR_BYTES),
        ("uniprot_sidecar", uniprot_sidecar_path, MAX_SIDECAR_BYTES),
    ):
        if spath is not None:
            p = Path(spath)
            sz = check_file_size_limit(
                p,
                max_bytes=max_b,
                label=f"Sidecar ({key})",
                empty_error_code="SIDECAR_MISMATCH",
                stage=stage,
            )
            total_bytes += sz
            if total_bytes > MAX_TOTAL_INPUT_BYTES:
                raise NBFoldError(
                    "UNSUPPORTED_INPUT_SIZE",
                    f"Total input size ({total_bytes} bytes) exceeds 5 MB limit ({MAX_TOTAL_INPUT_BYTES} bytes).",
                    stage=stage,
                )
            raw_b = p.read_bytes()
            sidecar_sha256s[key] = sha256_bytes(raw_b)
            try:
                parsed_json = json.loads(raw_b.decode("utf-8"))
            except Exception as exc:
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"Sidecar '{key}' is not valid UTF-8 JSON.",
                    stage=stage,
                ) from exc
            if key == "confidence_sidecar":
                confidence_sidecar_data = parsed_json
            elif key == "pae_sidecar":
                pae_sidecar_data = parsed_json
            elif key == "export_metadata":
                export_metadata_data = parsed_json
            elif key == "interpro_sidecar":
                if not isinstance(parsed_json, dict):
                    raise NBFoldError(
                        "SIDECAR_MISMATCH",
                        "InterPro sidecar must be a JSON object.",
                        stage=stage,
                    )
                interpro_sidecar_data = parsed_json
            elif key == "uniprot_sidecar":
                if not isinstance(parsed_json, dict):
                    raise NBFoldError(
                        "SIDECAR_MISMATCH",
                        "UniProt sidecar must be a JSON object.",
                        stage=stage,
                    )
                uniprot_sidecar_data = parsed_json

    if total_bytes > MAX_TOTAL_INPUT_BYTES:
        raise NBFoldError(
            "UNSUPPORTED_INPUT_SIZE",
            f"Total input size ({total_bytes} bytes) exceeds 5 MB limit ({MAX_TOTAL_INPUT_BYTES} bytes).",
            stage=stage,
        )

    # 2. Parse protein FASTA and structure format
    parsed_protein = parse_protein_fasta_file(Path(protein_fasta_path), stage=stage)
    fmt, structure_text, structure_sha256, _ = resolve_and_validate_structure_format(
        Path(structure_file_path), declared_format=declared_format, stage=stage
    )

    # 3. Pre-scan raw B-factor fields for NaN/Inf before Biopython parsing
    _prescan_nonfinite_bfactors(structure_text, fmt, stage=stage)

    # 4. Parse structure with Biopython
    try:
        parser = (
            PDBParser(PERMISSIVE=False, QUIET=True) if fmt == "pdb" else MMCIFParser(QUIET=True)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure("nbfold_intake", StringIO(structure_text))
    except Exception as exc:
        raise NBFoldError(
            "INVALID_STRUCTURE_FORMAT",
            "Structure file could not be parsed by strict PDB/mmCIF parser.",
            stage=stage,
        ) from exc

    models = list(structure.get_models())
    if model_number < 1 or model_number > len(models):
        raise NBFoldError(
            "MODEL_NOT_FOUND",
            f"Requested model {model_number} was not found (file contains {len(models)} model(s)).",
            stage=stage,
        )
    selected_model = models[model_number - 1]
    if chain_id not in selected_model:
        raise NBFoldError(
            "CHAIN_NOT_FOUND",
            f"Requested chain '{chain_id}' was not found in model {model_number}.",
            stage=stage,
        )
    chain = selected_model[chain_id]

    # 5. Inspect residues strictly for hetero/modified residues, disorder, sequence match, and CA atoms
    polymer_residues = []
    for res in chain.get_residues():
        hetflag = res.id[0]
        resname = res.resname.strip().upper()
        if hetflag == "W" or resname in {"HOH", "WAT"}:
            continue
        if hetflag != " ":
            raise NBFoldError(
                "MODIFIED_OR_AMBIGUOUS_RESIDUE",
                f"Chain '{chain_id}' contains hetero or modified residue '{resname}' (hetflag '{hetflag}'). "
                "Automatic repair or substitution is not supported.",
                stage=stage,
            )
        if res.is_disordered() != 0 or any(
            atom.is_disordered() != 0 for atom in res.get_atoms()
        ):
            raise NBFoldError(
                "DISORDERED_RESIDUE",
                f"Chain '{chain_id}' contains disordered residue or atom at {res.id}.",
                stage=stage,
            )
        if resname not in CANONICAL_3TO1:
            raise NBFoldError(
                "MODIFIED_OR_AMBIGUOUS_RESIDUE",
                f"Chain '{chain_id}' contains non-standard residue '{resname}' at {res.id}.",
                stage=stage,
            )
        polymer_residues.append(res)

    observed_sequence = "".join(
        CANONICAL_3TO1[r.resname.strip().upper()] for r in polymer_residues
    )
    if observed_sequence != parsed_protein.sequence:
        raise NBFoldError(
            "SEQUENCE_MISMATCH",
            f"Selected structure chain '{chain_id}' (length {len(observed_sequence)} aa) does not "
            f"exactly match input protein (length {len(parsed_protein.sequence)} aa). "
            "Exact full-chain residue identity and length match are required.",
            stage=stage,
            details={
                "observed_length_aa": len(observed_sequence),
                "expected_length_aa": len(parsed_protein.sequence),
                "chain_id": chain_id,
                "model_number": model_number,
            },
        )

    for idx, res in enumerate(polymer_residues, start=1):
        if "CA" not in res:
            raise NBFoldError(
                "MISSING_CA_ATOM",
                f"Residue {idx} ({res.resname.strip()}) in chain '{chain_id}' is missing a required CA atom.",
                stage=stage,
                details={"protein_position": idx, "chain_id": chain_id},
            )

    # 6. Resolve provider family and source basis
    resolved_provider, source_basis, export_version = _resolve_provider_and_source_basis(
        declared_provider=provider_family,
        structure_text=structure_text,
        export_metadata=export_metadata_data,
        stage=stage,
    )
    is_af = resolved_provider in {"af2", "af3"}

    # 7. Extract CA confidence values and validate finiteness/range
    ca_bfactors: list[float] = []
    atom_variance_flags: list[bool] = []
    for idx, res in enumerate(polymer_residues, start=1):
        atom_bfs = []
        for atom in res.get_atoms():
            val = float(atom.get_bfactor())
            if not math.isfinite(val):
                raise NBFoldError(
                    "INVALID_CONFIDENCE_VALUE",
                    f"Non-finite atom B-factor/confidence at residue {idx} in chain '{chain_id}'.",
                    stage=stage,
                )
            atom_bfs.append(val)
        ca_val = float(res["CA"].get_bfactor())
        if not math.isfinite(ca_val) or ca_val < 0.0 or (is_af and ca_val > 100.0):
            raise NBFoldError(
                "INVALID_CONFIDENCE_VALUE",
                f"CA confidence/B-factor ({ca_val}) at residue {idx} is outside valid bounds "
                f"({'[0, 100]' if is_af else '[0, inf)'}).",
                stage=stage,
                details={"protein_position": idx, "observed_value": ca_val},
            )
        ca_bfactors.append(ca_val)
        varies = (max(atom_bfs) - min(atom_bfs)) > 1e-3 if len(atom_bfs) > 1 else False
        atom_variance_flags.append(varies)

    # 8. Validate optional sidecars
    sidecar_obs, domain_rel, sidecar_token_indices, chain_pae_matrix = _validate_and_map_sidecars(
        confidence_sidecar_data=confidence_sidecar_data,
        pae_sidecar_data=pae_sidecar_data,
        chain_id=chain_id,
        protein_length=len(parsed_protein.sequence),
        ca_bfactors=ca_bfactors,
        stage=stage,
    )

    # 9. Reuse auto_annotator heuristic for candidate regions (unchanged objective/threshold)
    base_report = recommend_structure_pause_regions(
        protein=parsed_protein.sequence,
        structure_content=structure_text,
        format=fmt,
        chain_id=chain_id,
        model_number=model_number,
        is_alphafold=is_af,
        min_region_len=min_region_len,
        max_region_len=max_region_len,
    )
    recommended_positions = list(base_report.recommended_positions)

    # 10. Phase C structural annotation callers (domain_caller, tm_signal_caller, backbone_secondary_structure)
    domain_caller_report, domain_rel = analyze_domain_annotations(
        protein_sequence=parsed_protein.sequence,
        protein_sha256=parsed_protein.protein_sha256,
        ca_bfactors=ca_bfactors,
        pae_matrix=chain_pae_matrix,
        interpro_sidecar_data=interpro_sidecar_data,
        interpro_sidecar_sha256=sidecar_sha256s["interpro_sidecar"],
        recommended_positions=recommended_positions,
        base_domain_relative_confidence=domain_rel,
        stage=stage,
    )
    tm_signal_caller_report = analyze_uniprot_tm_signal_annotations(
        protein_sequence=parsed_protein.sequence,
        protein_sha256=parsed_protein.protein_sha256,
        uniprot_sidecar_data=uniprot_sidecar_data,
        uniprot_sidecar_sha256=sidecar_sha256s["uniprot_sidecar"],
        recommended_positions=recommended_positions,
        stage=stage,
    )
    backbone_ss_report = analyze_backbone_secondary_structure(
        polymer_residues=polymer_residues,
        protein_sequence=parsed_protein.sequence,
        recommended_positions=recommended_positions,
    )

    candidate_regions_list = [
        {
            "start_1based": r.start_1based,
            "end_1based": r.end_1based,
            "length": r.length,
            "region_type": r.region_type,
            "reason": r.reason,
            "heuristic_basis": r.heuristic_basis,
            "hypothesis_class": "unconfirmed_structural_observation_hypothesis",
        }
        for r in base_report.candidate_regions
    ]
    structural_review_summary = build_structural_review_summary(
        protein_length=len(parsed_protein.sequence),
        candidate_regions=candidate_regions_list,
        recommended_positions=recommended_positions,
        domain_caller_report=domain_caller_report,
        tm_signal_caller_report=tm_signal_caller_report,
        backbone_ss_report=backbone_ss_report,
    )

    if resolved_provider == "af3":
        confidence_metric = "CA_atom_confidence"
        confidence_basis = (
            "AlphaFold 3 atom-level confidence sampled at CA atom (0-100); "
            "not averaged across residue atoms and distinct from AF2 residue-level pLDDT"
        )
    elif resolved_provider == "af2":
        confidence_metric = "residue_plddt"
        confidence_basis = "AlphaFold 2 per-residue pLDDT recorded on CA atom (0-100)"
    else:
        confidence_metric = "experimental_bfactor"
        confidence_basis = "Experimental crystallographic/NMR B-factors on CA atom (not calibrated pLDDT)"

    mmcif_label_map = (
        _extract_mmcif_label_seq_ids(structure_text, chain_id) if fmt == "mmcif" else {}
    )
    conf_source_label = (
        "structure_bfactor_ca+verified_sidecar"
        if sidecar_obs.get("confidence_sidecar_verified")
        else "structure_bfactor_ca"
    )

    position_mapping: list[dict[str, Any]] = []
    for idx, (res, aa, ca_bf, varies, tok_idx) in enumerate(
        zip(
            polymer_residues,
            parsed_protein.sequence,
            ca_bfactors,
            atom_variance_flags,
            sidecar_token_indices,
        ),
        start=1,
    ):
        auth_seq = int(res.id[1])
        ins_code = res.id[2].strip()
        lbl_seq = mmcif_label_map.get(auth_seq, idx if not ins_code else None)
        position_mapping.append(
            {
                "protein_position": idx,
                "amino_acid": aa,
                "resname_3letter": res.resname.strip().upper(),
                "chain_id": chain_id,
                "model_number": model_number,
                "auth_seq_id": auth_seq,
                "insertion_code": ins_code,
                "label_seq_id": lbl_seq,
                "selected_atom": "CA",
                "confidence_value": round(ca_bf, 2),
                "confidence_metric": confidence_metric,
                "confidence_source": conf_source_label,
                "atom_confidence_varies_within_residue": varies,
                "sidecar_token_index": tok_idx,
            }
        )

    status = base_report.status  # "AVAILABLE" or "NO_CANDIDATES"
    workflow_state = (
        "AWAITING_USER_REVIEW" if status == "AVAILABLE" else "NEEDS_ANNOTATION"
    )

    active_sidecar_sha256s: dict[str, str | None] = {
        "confidence_sidecar": sidecar_sha256s["confidence_sidecar"],
        "pae_sidecar": sidecar_sha256s["pae_sidecar"],
    }
    if sidecar_sha256s["interpro_sidecar"] is not None:
        active_sidecar_sha256s["interpro_sidecar"] = sidecar_sha256s["interpro_sidecar"]
    if sidecar_sha256s["uniprot_sidecar"] is not None:
        active_sidecar_sha256s["uniprot_sidecar"] = sidecar_sha256s["uniprot_sidecar"]

    annotation_payload: dict[str, Any] = {
        "schema_version": ANNOTATION_SCHEMA_VERSION,
        "pipeline_version": NBFOLD_PIPELINE_VERSION,
        "recommender_version": RECOMMENDER_VERSION,
        "status": status,
        "workflow_state": workflow_state,
        "reason": base_report.reason,
        "reviewed": False,  # Hard invariant: always False
        "protein_length": len(parsed_protein.sequence),
        "protein_sha256": parsed_protein.protein_sha256,
        "protein_fasta_file_sha256": parsed_protein.file_sha256,
        "structure_sha256": structure_sha256,
        "sidecar_sha256s": active_sidecar_sha256s,
        "export_metadata_sha256": sidecar_sha256s["export_metadata"],
        "format": fmt,
        "chain_id": chain_id,
        "model_number": model_number,
        "exact_sequence_match": True,
        "provider_family": resolved_provider,
        "source_basis": source_basis,
        "export_version": export_version,
        "provenance_verification": "unverified_local_claim",
        "is_alphafold_provenance": is_af,
        "confidence_metric": confidence_metric,
        "confidence_basis": confidence_basis,
        "selected_atom_policy": "CA_atom_only",
        "heuristic_config": {
            "plddt_low_confidence_threshold": 70.0 if is_af else None,
            "experimental_relative_bfactor_multiplier": 1.2 if not is_af else None,
            "experimental_min_bfactor": 20.0 if not is_af else None,
            "min_region_len": min_region_len,
            "max_region_len": max_region_len,
            "segmentation_mode": "contiguous_bounded_windows",
            "no_midpoint_fallback": True,
        },
        "candidate_regions": candidate_regions_list,
        "recommended_positions": recommended_positions,
        "coverage": base_report.coverage,
        "position_mapping": position_mapping,
        "sidecar_observations": sidecar_obs,
        "domain_relative_confidence": domain_rel,
        "independent_structural_tools": {
            "domain_boundary_assignment": domain_caller_report,
            "transmembrane_topology": tm_signal_caller_report,
            "backbone_secondary_structure": backbone_ss_report,
            "secondary_structure_dssp": {
                "status": "NOT_AVAILABLE",
                "tool": None,
                "reason": (
                    "External licensed DSSP binary is not executed; see backbone_secondary_structure "
                    "for local N/CA/C/O dihedral and H-bond geometry observations."
                ),
            },
        },
        "structural_review_summary": structural_review_summary,
        "default_config_digest": default_config_digest,
        "resolved_config_preview": resolved_config_preview,
        "limitations": LIMITATIONS_STATEMENT,
    }

    content_hash = compute_annotation_content_sha256(annotation_payload)
    annotation_payload["annotation_content_sha256"] = content_hash
    annotation_payload["review_manifest_template"] = {
        "schema_version": "nbfold-review-v1",
        "review_type": "LOCAL_USER_ACKNOWLEDGEMENT",
        "positions_reviewed": False,
        "reviewed_at": None,
        "reviewer_label": None,
        "protein_sha256": parsed_protein.protein_sha256,
        "structure_sha256": structure_sha256,
        "sidecar_sha256s": dict(active_sidecar_sha256s),
        "chain_id": chain_id,
        "model_number": model_number,
        "annotation_sha256": content_hash,
        "reviewed_positions": list(recommended_positions),
        "config_digest": default_config_digest,
        "conflicts_acknowledged": not structural_review_summary["has_conflicts"],
        "conflict_resolution_rationale": None,
    }

    return annotation_payload
