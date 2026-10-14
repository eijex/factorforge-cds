"""Phase C structural annotation callers and conflict analysis for NB-Fold (Job 399).

Implements PRD-JOB-399-NB-FOLD-PHASE-C-v1.0.0:
1. `domain_caller`: PAE sliding-window domain segmentation (adapted from Science
   `alphafold-database-fetch-and-analyze` `analyze_pae.py`) and/or local InterPro
   domain sidecar integration, computing inter-domain linkers, terminal tails,
   per-domain `domain_relative_confidence`, and boundary/core conflicts.
2. `tm_signal_caller`: Local UniProt feature sidecar integration for Signal peptide,
   Transmembrane, Intramembrane, Topological domain, Disulfide bond, and Glycosylation
   annotations with strict sequence/length/coordinate validation and conflict detection.
3. `backbone_secondary_structure`: Local backbone (`N`, `CA`, `C`, `O`) dihedral
   (phi, psi) and H-bond/CA geometry observer (`H`, `E`, `C`), strictly separated from
   external `dssp_secondary_structure` (`NOT_AVAILABLE`) and reporting `NOT_AVAILABLE`
   on CA-only structures without inventing secondary structure.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from factorforge.nbfold.contracts import (
    BACKBONE_SS_CALLER_VERSION,
    DOMAIN_CALLER_VERSION,
    TM_SIGNAL_CALLER_VERSION,
    NBFoldError,
    sha256_text,
)

MIN_PAE_SEGMENTATION_PROTEIN_AA = 20
DEFAULT_PAE_DISTANCE_CUTOFF_A = 7.0
DEFAULT_PAE_MERGE_CUTOFF_A = 12.0
DEFAULT_PAE_MIN_SUBDOMAIN_AA = 8
DEFAULT_PAE_MIN_DOMAIN_AA = 8

PAE_PARAMETER_PROFILES: dict[str, dict[str, Any]] = {
    "v1_2_1_baseline": {
        "profile_id": "v1_2_1_baseline",
        "description": (
            "Job 399 v1.2.1 baseline PAE segmentation parameters "
            "(preserved as default until public benchmark comparison review)."
        ),
        "min_protein_aa": MIN_PAE_SEGMENTATION_PROTEIN_AA,
        "distance_cutoff_angstrom": DEFAULT_PAE_DISTANCE_CUTOFF_A,
        "merge_cutoff_angstrom": DEFAULT_PAE_MERGE_CUTOFF_A,
        "min_subdomain_aa": DEFAULT_PAE_MIN_SUBDOMAIN_AA,
        "min_domain_aa": DEFAULT_PAE_MIN_DOMAIN_AA,
    },
    "prd_v1_0_0_candidate": {
        "profile_id": "prd_v1_0_0_candidate",
        "description": (
            "Original PRD-JOB-399 v1.0.0 §4.1 candidate parameters "
            "(Science analyze_pae.py defaults with post-merge min_domain_aa=40)."
        ),
        "min_protein_aa": MIN_PAE_SEGMENTATION_PROTEIN_AA,
        "distance_cutoff_angstrom": 7.0,
        "merge_cutoff_angstrom": 15.0,
        "min_subdomain_aa": 10,
        "min_domain_aa": 40,
    },
}

_INTERPRO_TYPE_ALIASES: dict[str, str] = {
    "domain": "domain",
    "family": "family",
    "homologous_superfamily": "homologous_superfamily",
    "superfamily": "homologous_superfamily",
    "repeat": "repeat",
    "site": "site",
    "conserved_site": "site",
    "active_site": "site",
    "binding_site": "site",
    "ptm": "site",
    "unclassified": "unclassified",
}



def _validate_sidecar_binding(data: Mapping[str, Any], protein: str, digest: str, *, label: str, stage: str) -> None:
    """Validate supplied binding fields without coercing malformed lengths or identities."""
    found = False
    for field in ("protein_sha256", "sequence", "protein_length"):
        if field not in data:
            continue
        value = data[field]
        if field == "sequence" and isinstance(value, Mapping):
            fields = [("sequence.value", value["value"])] if "value" in value else []
            if "length" in value:
                fields.append(("protein_length", value["length"]))
        else:
            fields = [(field, value)]
        for name, item in fields:
            valid = False
            if name == "protein_sha256":
                valid = isinstance(item, str) and item == digest
            elif name == "protein_length":
                valid = type(item) is int and item == len(protein)
            elif name in {"sequence", "sequence.value"}:
                valid = isinstance(item, str) and "".join(item.upper().split()) == protein
            if not valid:
                raise NBFoldError("SIDECAR_MISMATCH", f"{label} sidecar has invalid or mismatched {name}.", stage=stage)
            found = True
    if not found:
        raise NBFoldError("SIDECAR_MISMATCH", f"{label} sidecar requires sequence, protein_sha256, or protein_length binding.", stage=stage)


def _vec_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _vec_dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vec_norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(_vec_dot(a, a))


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return _vec_norm(_vec_sub(a, b))


def _dihedral_degrees(
    p0: tuple[float, float, float],
    p1: tuple[float, float, float],
    p2: tuple[float, float, float],
    p3: tuple[float, float, float],
) -> float | None:
    """Compute IUPAC dihedral angle in degrees [-180, 180] for four 3D points."""
    b0 = _vec_sub(p0, p1)
    b1 = _vec_sub(p2, p1)
    b2 = _vec_sub(p3, p2)
    b1_norm = _vec_norm(b1)
    if b1_norm < 1e-6:
        return None
    b1_unit = (b1[0] / b1_norm, b1[1] / b1_norm, b1[2] / b1_norm)
    v = (
        b0[0] - _vec_dot(b0, b1_unit) * b1_unit[0],
        b0[1] - _vec_dot(b0, b1_unit) * b1_unit[1],
        b0[2] - _vec_dot(b0, b1_unit) * b1_unit[2],
    )
    w = (
        b2[0] - _vec_dot(b2, b1_unit) * b1_unit[0],
        b2[1] - _vec_dot(b2, b1_unit) * b1_unit[1],
        b2[2] - _vec_dot(b2, b1_unit) * b1_unit[2],
    )
    if _vec_norm(v) < 1e-6 or _vec_norm(w) < 1e-6:
        return None
    x = _vec_dot(v, w)
    y = _vec_dot(_vec_cross(b1_unit, v), w)
    deg = math.degrees(math.atan2(y, x))
    return round(deg, 2) if math.isfinite(deg) else None


def _extract_atom_coord(res: Any, atom_name: str) -> tuple[float, float, float] | None:
    if atom_name not in res:
        return None
    coord = res[atom_name].get_coord()
    if len(coord) != 3:
        return None
    x, y, z = float(coord[0]), float(coord[1]), float(coord[2])
    if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
        return None
    return (x, y, z)


def analyze_backbone_secondary_structure(
    *,
    polymer_residues: Sequence[Any],
    protein_sequence: str,
    recommended_positions: Sequence[int],
) -> dict[str, Any]:
    """Observe local secondary structure (`H`, `E`, `C`) from full `N, CA, C, O` backbone coordinates.

    Returns `status="NOT_AVAILABLE"` when no residue has complete `N, CA, C, O` backbone atoms
    (e.g., CA-only or CA/CB-only fixtures), never fabricating secondary structure.
    """
    n_res = len(polymer_residues)
    backbone_coords: list[dict[str, tuple[float, float, float] | None]] = []
    complete_count = 0

    for res in polymer_residues:
        atoms = {
            name: _extract_atom_coord(res, name)
            for name in ("N", "CA", "C", "O")
        }
        if all(atoms[k] is not None for k in ("N", "CA", "C", "O")):
            complete_count += 1
        backbone_coords.append(atoms)

    if complete_count == 0:
        return {
            "status": "NOT_AVAILABLE",
            "source": "structure_backbone_coordinates",
            "method": "local_backbone_dihedral_hbond_v1",
            "version": BACKBONE_SS_CALLER_VERSION,
            "reason": (
                "Structure lacks full N, CA, C, O backbone coordinates (CA/CB-only trace); "
                "secondary structure is not inferred without full backbone geometry."
            ),
            "complete_backbone_residues": 0,
            "total_residues": n_res,
            "residue_assignments": [],
            "segments": [],
            "conflicts": [],
        }

    status = "AVAILABLE" if complete_count == n_res else "PARTIAL"

    # Compute phi, psi and initial raw classification per residue
    phis: list[float | None] = [None] * n_res
    psis: list[float | None] = [None] * n_res
    raw_states: list[str] = ["U"] * n_res
    uncertainties: list[str] = ["missing_backbone_atoms"] * n_res

    for i in range(n_res):
        cur = backbone_coords[i]
        if any(cur[k] is None for k in ("N", "CA", "C", "O")):
            raw_states[i] = "U"
            uncertainties[i] = "missing_backbone_atoms"
            continue

        # phi_i = dihedral(C_{i-1}, N_i, CA_i, C_i)
        if i > 0 and backbone_coords[i - 1]["C"] is not None:
            c_prev = backbone_coords[i - 1]["C"]
            n_cur = cur["N"]
            if c_prev is not None and n_cur is not None and _dist(c_prev, n_cur) <= 2.2:
                phis[i] = _dihedral_degrees(c_prev, n_cur, cur["CA"], cur["C"])  # type: ignore[arg-type]

        # psi_i = dihedral(N_i, CA_i, C_i, N_{i+1})
        if i + 1 < n_res and backbone_coords[i + 1]["N"] is not None:
            c_cur = cur["C"]
            n_next = backbone_coords[i + 1]["N"]
            if c_cur is not None and n_next is not None and _dist(c_cur, n_next) <= 2.2:
                psis[i] = _dihedral_degrees(cur["N"], cur["CA"], c_cur, n_next)  # type: ignore[arg-type]

        phi = phis[i]
        psi = psis[i]
        if phi is None and psi is None:
            raw_states[i] = "U"
            uncertainties[i] = "unavailable_backbone_dihedrals"
            continue

        # Check i -> i+4 / i+3 carbonyl-amide H-bond distance or CA(i)-CA(i+3) helical pitch
        has_helical_hbond = False
        for offset in (3, 4):
            if i + offset < n_res:
                o_i = cur["O"]
                n_j = backbone_coords[i + offset]["N"]
                if o_i is not None and n_j is not None and _dist(o_i, n_j) <= 3.6:
                    has_helical_hbond = True
            if i - offset >= 0:
                o_prev = backbone_coords[i - offset]["O"]
                n_i = cur["N"]
                if o_prev is not None and n_i is not None and _dist(o_prev, n_i) <= 3.6:
                    has_helical_hbond = True

        if phi is not None and psi is not None:
            in_helix_basin = (-140.0 <= phi <= -30.0) and (-80.0 <= psi <= 20.0)
            in_core_helix = (-100.0 <= phi <= -45.0) and (-70.0 <= psi <= -20.0)
            in_strand_basin = (phi <= -70.0 or phi >= 150.0) and (psi >= 80.0 or psi <= -140.0)
            in_core_strand = (-165.0 <= phi <= -90.0) and (95.0 <= psi <= 170.0)

            if in_helix_basin:
                raw_states[i] = "H"
                uncertainties[i] = (
                    "high_confidence_geometry"
                    if (in_core_helix or has_helical_hbond)
                    else "boundary_transition"
                )
            elif in_strand_basin:
                raw_states[i] = "E"
                uncertainties[i] = (
                    "high_confidence_geometry" if in_core_strand else "boundary_transition"
                )
            else:
                raw_states[i] = "C"
                uncertainties[i] = "high_confidence_geometry"
        else:
            # Terminal residue with only one dihedral angle
            if has_helical_hbond and (
                (psi is not None and -80.0 <= psi <= 20.0)
                or (phi is not None and -140.0 <= phi <= -30.0)
            ):
                raw_states[i] = "H"
            elif (psi is not None and (psi >= 90.0 or psi <= -140.0)) or (
                phi is not None and (phi <= -85.0 or phi >= 150.0)
            ):
                raw_states[i] = "E"
            else:
                raw_states[i] = "C"
            uncertainties[i] = "terminal_residue_single_dihedral"

    # Smooth isolated 1-residue/2-residue H (<3) and isolated 1-residue E (<2) into C
    smoothed = list(raw_states)
    idx = 0
    while idx < n_res:
        st = smoothed[idx]
        j = idx
        while j < n_res and smoothed[j] == st:
            j += 1
        run_len = j - idx
        if (st == "H" and run_len < 3) or (st == "E" and run_len < 2):
            for k in range(idx, j):
                smoothed[k] = "C"
                uncertainties[k] = "isolated_dihedral_smoothed_to_coil"
        idx = j

    # Build contiguous segments and residue_assignments
    segments: list[dict[str, Any]] = []
    idx = 0
    while idx < n_res:
        st = smoothed[idx]
        j = idx
        while j < n_res and smoothed[j] == st:
            j += 1
        segments.append(
            {
                "ss_state": st,
                "start_1based": idx + 1,
                "end_1based": j,
                "length": j - idx,
            }
        )
        # Mark immediate segment boundary residues for H/E
        if st in {"H", "E"} and (j - idx) >= 2:
            if uncertainties[idx] == "high_confidence_geometry":
                uncertainties[idx] = "boundary_transition"
            if uncertainties[j - 1] == "high_confidence_geometry":
                uncertainties[j - 1] = "boundary_transition"
        idx = j

    residue_assignments: list[dict[str, Any]] = []
    for i in range(n_res):
        residue_assignments.append(
            {
                "protein_position": i + 1,
                "amino_acid": protein_sequence[i],
                "ss_state": smoothed[i],
                "phi_deg": phis[i],
                "psi_deg": psis[i],
                "uncertainty": uncertainties[i],
            }
        )

    assigned_count = sum(state != "U" for state in smoothed)
    status = "NOT_AVAILABLE" if assigned_count == 0 else ("AVAILABLE" if assigned_count == n_res else "PARTIAL")

    conflicts = detect_secondary_structure_conflicts(
        segments=segments,
        positions=recommended_positions,
    )

    return {
        "status": status,
        "source": "structure_backbone_coordinates",
        "method": "local_backbone_dihedral_hbond_v1",
        "version": BACKBONE_SS_CALLER_VERSION,
        "reason": None,
        "complete_backbone_residues": complete_count,
        "total_residues": n_res,
        "residue_assignments": residue_assignments,
        "segments": segments,
        "conflicts": conflicts,
    }


def detect_secondary_structure_conflicts(
    *,
    segments: Sequence[Mapping[str, Any]],
    positions: Sequence[int],
) -> list[dict[str, Any]]:
    """Flag positions that fall inside the interior core of a rigid Helix (len>=4) or Strand (len>=3)."""
    pos_set = {int(p) for p in positions}
    conflicts: list[dict[str, Any]] = []
    for seg in segments:
        st = seg.get("ss_state")
        start = int(seg["start_1based"])
        end = int(seg["end_1based"])
        length = int(seg["length"])
        if (st == "H" and length >= 4) or (st == "E" and length >= 3):
            core_positions = [p for p in range(start + 1, end) if p in pos_set]
            if core_positions:
                ss_label = "alpha-helix (H)" if st == "H" else "beta-strand (E)"
                conflicts.append(
                    {
                        "conflict_code": "CANDIDATE_IN_RIGID_SECONDARY_STRUCTURE_CORE",
                        "source": "backbone_secondary_structure",
                        "severity": "WARNING",
                        "ss_state": st,
                        "segment_start_1based": start,
                        "segment_end_1based": end,
                        "overlapping_positions": sorted(core_positions),
                        "message": (
                            f"Selected position(s) {sorted(core_positions)} fall inside the interior core "
                            f"of a {ss_label} segment ({start}-{end})."
                        ),
                    }
                )
    return conflicts


def _segment_pae_domains(
    pae_matrix: Sequence[Sequence[float]],
    ca_bfactors: Sequence[float],
    *,
    distance_cutoff: float = DEFAULT_PAE_DISTANCE_CUTOFF_A,
    merge_cutoff: float = DEFAULT_PAE_MERGE_CUTOFF_A,
    min_subdomain_aa: int = DEFAULT_PAE_MIN_SUBDOMAIN_AA,
    min_domain_aa: int = DEFAULT_PAE_MIN_DOMAIN_AA,
    min_protein_aa: int = MIN_PAE_SEGMENTATION_PROTEIN_AA,
) -> list[dict[str, Any]]:
    """Sliding-window PAE domain segmentation adapted from Science `analyze_pae.py`."""
    n_res = len(pae_matrix)
    if n_res < min_protein_aa:
        return []

    subdomains: list[list[int]] = []
    current: list[int] = []
    for i in range(n_res):
        if not current:
            current.append(i)
            continue
        win = min(20, len(current))
        recent = current[-win:]
        pae_sum = sum(float(pae_matrix[r][i]) + float(pae_matrix[i][r]) for r in recent)
        avg_pae = pae_sum / (2.0 * win)
        if avg_pae < distance_cutoff:
            current.append(i)
        else:
            if len(current) >= min_subdomain_aa:
                subdomains.append(current)
            current = [i]
    if len(current) >= min_subdomain_aa:
        subdomains.append(current)

    raw_bounds = [[comp[0] + 1, comp[-1] + 1] for comp in subdomains]
    if not raw_bounds:
        return []

    merged: list[list[int]] = [raw_bounds[0]]
    for b in raw_bounds[1:]:
        prev_start, prev_end = merged[-1][0] - 1, merged[-1][1] - 1
        cur_start, cur_end = b[0] - 1, b[1] - 1
        lookback = max(prev_start, prev_end - 30)
        lookfwd = min(cur_end, cur_start + 30)
        pae_sum = 0.0
        n_pairs = 0
        for r1 in range(lookback, prev_end + 1):
            for r2 in range(cur_start, lookfwd + 1):
                pae_sum += float(pae_matrix[r1][r2]) + float(pae_matrix[r2][r1])
                n_pairs += 2
        if n_pairs > 0 and (pae_sum / n_pairs) < merge_cutoff:
            merged[-1][1] = b[1]
        else:
            merged.append(b)

    filtered_merged = [
        (s, e) for s, e in merged if (e - s + 1) >= min_domain_aa
    ]

    domains: list[dict[str, Any]] = []
    for idx, (s, e) in enumerate(filtered_merged, start=1):
        vals = [
            float(pae_matrix[r1][r2])
            for r1 in range(s - 1, e)
            for r2 in range(s - 1, e)
        ]
        bfs = [float(ca_bfactors[r]) for r in range(s - 1, e)]
        mean_intra = round(sum(vals) / len(vals), 3) if vals else 0.0
        max_intra = round(max(vals), 3) if vals else 0.0
        domains.append(
            {
                "domain_id": f"PAE_D{idx}",
                "source": "pae_sliding_window",
                "version": DOMAIN_CALLER_VERSION,
                "start_1based": s,
                "end_1based": e,
                "length": e - s + 1,
                "mean_intra_pae_angstrom": mean_intra,
                "max_intra_pae_angstrom": max_intra,
                "mean_ca_confidence": round(sum(bfs) / len(bfs), 2) if bfs else 0.0,
                "min_ca_confidence": round(min(bfs), 2) if bfs else 0.0,
                "max_ca_confidence": round(max(bfs), 2) if bfs else 0.0,
                "uncertainty": {
                    "metric": "mean_intra_pae_angstrom",
                    "value": mean_intra,
                    "distance_cutoff_angstrom": distance_cutoff,
                    "merge_cutoff_angstrom": merge_cutoff,
                    "min_protein_aa": min_protein_aa,
                    "min_subdomain_aa": min_subdomain_aa,
                    "min_domain_aa": min_domain_aa,
                },
            }
        )
    return domains


def _parse_interpro_sidecar(
    interpro_data: Mapping[str, Any],
    *,
    protein_sequence: str,
    protein_sha256: str,
    interpro_sidecar_sha256: str | None,
    stage: str,
) -> list[dict[str, Any]]:
    """Validate and normalize local InterPro sidecar JSON, preserving entry_type."""
    n_res = len(protein_sequence)
    _validate_sidecar_binding(interpro_data, protein_sequence, protein_sha256, label="InterPro", stage=stage)
    if (
        "protein_sha256" in interpro_data
        and interpro_data["protein_sha256"] is not None
        and str(interpro_data["protein_sha256"]) != protein_sha256
    ):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            "InterPro sidecar protein_sha256 does not match input protein sequence.",
            stage=stage,
        )
    if "sequence" in interpro_data and isinstance(interpro_data["sequence"], str):
        clean_seq = "".join(interpro_data["sequence"].upper().split())
        if clean_seq != protein_sequence:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "InterPro sidecar sequence does not match input protein sequence.",
                stage=stage,
            )
    if (
        "protein_length" in interpro_data
        and interpro_data["protein_length"] is not None
        and int(interpro_data["protein_length"]) != n_res
    ):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            f"InterPro sidecar protein_length ({interpro_data['protein_length']}) does not match input protein ({n_res} aa).",
            stage=stage,
        )

    raw_entries = interpro_data.get("entries", interpro_data.get("matches"))
    if raw_entries is None and isinstance(interpro_data.get("results"), list):
        raw_entries = interpro_data["results"]
    if raw_entries is None and isinstance(interpro_data.get("domains"), list):
        raw_entries = interpro_data["domains"]
    if not isinstance(raw_entries, list):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            "InterPro sidecar must contain an 'entries', 'results', or 'domains' list.",
            stage=stage,
        )

    annotations: list[dict[str, Any]] = []
    for entry in raw_entries:
        if not isinstance(entry, Mapping):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "InterPro sidecar entry must be a JSON object.",
                stage=stage,
            )
        meta = entry.get("metadata") if isinstance(entry.get("metadata"), Mapping) else entry
        accession = str(meta.get("accession") or entry.get("accession") or "").strip()
        source_db = str(
            meta.get("source_database")
            or entry.get("source_db")
            or entry.get("database")
            or "interpro"
        ).strip().lower()
        name = str(meta.get("name") or entry.get("name") or accession).strip()
        version = str(
            entry.get("version")
            or meta.get("version")
            or interpro_data.get("interpro_version")
            or interpro_data.get("source_version")
            or "local_sidecar_v1"
        )
        if not accession:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "InterPro sidecar entry is missing required 'accession'.",
                stage=stage,
            )

        raw_type: Any = None
        for key in ("entry_type", "type"):
            if key in entry and entry[key] is not None:
                raw_type = entry[key]
                break
            if isinstance(entry.get("metadata"), Mapping) and key in entry["metadata"] and entry["metadata"][key] is not None:
                raw_type = entry["metadata"][key]
                break

        if raw_type is None:
            entry_type = "domain"
            entry_type_source = "default_untyped_entry"
        else:
            if not isinstance(raw_type, str):
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"InterPro entry '{accession}' type must be a string.",
                    stage=stage,
                )
            norm_type = raw_type.strip().lower().replace(" ", "_").replace("-", "_")
            if norm_type not in _INTERPRO_TYPE_ALIASES:
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"InterPro entry '{accession}' has unsupported entry type '{raw_type}'.",
                    stage=stage,
                )
            entry_type = _INTERPRO_TYPE_ALIASES[norm_type]
            entry_type_source = "explicit_entry_type"

        used_for_domain_consensus = entry_type == "domain"

        locations = entry.get("locations")
        if locations is None and ("start" in entry or "start_1based" in entry):
            locations = [
                {
                    "start": entry.get("start_1based", entry.get("start")),
                    "end": entry.get("end_1based", entry.get("end")),
                    "score": entry.get("score"),
                    "evalue": entry.get("evalue"),
                }
            ]
        if locations is None and isinstance(entry.get("proteins"), list) and entry["proteins"]:
            first_prot = entry["proteins"][0]
            if isinstance(first_prot, Mapping):
                prot_len = first_prot.get("protein_length")
                if prot_len is not None and (type(prot_len) is not int or prot_len != n_res):
                    raise NBFoldError(
                        "SIDECAR_MISMATCH",
                        f"InterPro entry protein_length ({prot_len}) does not match input protein ({n_res}).",
                        stage=stage,
                    )
                locations = []
                for loc in first_prot.get("entry_protein_locations", []) or []:
                    if isinstance(loc, Mapping):
                        for frag in loc.get("fragments", []) or []:
                            if isinstance(frag, Mapping):
                                locations.append(
                                    {
                                        "start": frag.get("start"),
                                        "end": frag.get("end"),
                                        "score": loc.get("score"),
                                        "evalue": loc.get("evalue"),
                                    }
                                )
        if not isinstance(locations, list) or not locations:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                f"InterPro entry '{accession}' has no valid residue locations.",
                stage=stage,
            )

        for loc_idx, loc in enumerate(locations, start=1):
            if not isinstance(loc, Mapping):
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"InterPro entry '{accession}' location must be an object.",
                    stage=stage,
                )
            start = loc.get("start")
            end = loc.get("end")
            if (
                type(start) is not int
                or type(end) is not int
                or not (1 <= start <= end <= n_res)
            ):
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    f"InterPro entry '{accession}' location [{start}, {end}] is outside 1..{n_res}.",
                    stage=stage,
                )
            evalue = loc.get("evalue")
            score = loc.get("score")
            for metric_name, metric in (("evalue", evalue), ("score", score)):
                if metric is not None and (
                    isinstance(metric, bool) or not isinstance(metric, (int, float))
                    or not math.isfinite(float(metric))
                ):
                    raise NBFoldError("SIDECAR_MISMATCH", f"InterPro {metric_name} must be finite numeric data or null.", stage=stage)

            annotations.append(
                {
                    "domain_id": f"{accession}_{loc_idx}",
                    "source": f"interpro:{source_db}",
                    "accession": accession,
                    "name": name,
                    "entry_type": entry_type,
                    "entry_type_source": entry_type_source,
                    "used_for_domain_consensus": used_for_domain_consensus,
                    "version": version,
                    "sidecar_sha256": interpro_sidecar_sha256,
                    "start_1based": start,
                    "end_1based": end,
                    "length": end - start + 1,
                    "evalue": float(evalue) if isinstance(evalue, (int, float)) else None,
                    "score": float(score) if isinstance(score, (int, float)) else None,
                    "uncertainty": {
                        "evalue": float(evalue) if isinstance(evalue, (int, float)) else None,
                        "score": float(score) if isinstance(score, (int, float)) else None,
                        "uncertainty_class": (
                            "curated_member_signature"
                            if evalue is not None or score is not None
                            else "unscored_interval_annotation"
                        ),
                    },
                }
            )
    annotations.sort(key=lambda d: (d["start_1based"], d["end_1based"], d["accession"]))
    return annotations


def analyze_domain_annotations(
    *,
    protein_sequence: str,
    protein_sha256: str,
    ca_bfactors: Sequence[float],
    pae_matrix: Sequence[Sequence[float]] | None,
    interpro_sidecar_data: Mapping[str, Any] | None,
    interpro_sidecar_sha256: str | None,
    recommended_positions: Sequence[int],
    base_domain_relative_confidence: Mapping[str, Any],
    stage: str = "annotate",
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Run Phase C `domain_caller` (PAE sliding-window + InterPro sidecar) and update `domain_relative_confidence`."""
    n_res = len(protein_sequence)
    pae_domains = (
        _segment_pae_domains(pae_matrix, ca_bfactors)
        if pae_matrix is not None
        else []
    )
    interpro_annotations = (
        _parse_interpro_sidecar(
            interpro_sidecar_data,
            protein_sequence=protein_sequence,
            protein_sha256=protein_sha256,
            interpro_sidecar_sha256=interpro_sidecar_sha256,
            stage=stage,
        )
        if interpro_sidecar_data is not None
        else []
    )
    interpro_domains = [
        dict(item) for item in interpro_annotations if item.get("used_for_domain_consensus", True)
    ]
    interpro_non_domain_annotations = [
        dict(item) for item in interpro_annotations if not item.get("used_for_domain_consensus", True)
    ]

    sources_used: list[str] = []
    if pae_domains:
        sources_used.append("pae_sliding_window")
    if interpro_annotations:
        sources_used.append("interpro_sidecar")

    domain_rel = dict(base_domain_relative_confidence)
    if not sources_used:
        reason = (
            "No InterPro domain sidecar supplied and PAE matrix absent or below minimum domain segmentation length "
            f"({MIN_PAE_SEGMENTATION_PROTEIN_AA} aa)."
        )
        return (
            {
                "status": "NOT_AVAILABLE",
                "tool": None,
                "version": DOMAIN_CALLER_VERSION,
                "sources_used": [],
                "reason": reason,
                "pae_parameter_profile": dict(PAE_PARAMETER_PROFILES["v1_2_1_baseline"]),
                "pae_domains": [],
                "interpro_domains": [],
                "interpro_non_domain_annotations": [],
                "interpro_annotations": [],
                "consensus_domains": [],
                "inter_domain_linkers": [],
                "terminal_regions": [],
                "conflicts": [],
            },
            domain_rel,
        )

    # Check cross-source conflicts if both PAE domains and InterPro domain-level entries exist
    conflicts: list[dict[str, Any]] = []
    if pae_domains and interpro_domains:
        # Compare primary non-overlapping intervals
        if len(pae_domains) != len(interpro_domains):
            conflicts.append(
                {
                    "conflict_code": "PAE_INTERPRO_DOMAIN_COUNT_MISMATCH",
                    "source": "domain_caller",
                    "severity": "WARNING",
                    "pae_domain_count": len(pae_domains),
                    "interpro_domain_count": len(interpro_domains),
                    "overlapping_positions": [],
                    "message": (
                        f"PAE segmentation detected {len(pae_domains)} domain(s) while "
                        f"InterPro sidecar reported {len(interpro_domains)} domain entry location(s)."
                    ),
                }
            )
        else:
            for p_dom, i_dom in zip(pae_domains, interpro_domains):
                start_diff = abs(int(p_dom["start_1based"]) - int(i_dom["start_1based"]))
                end_diff = abs(int(p_dom["end_1based"]) - int(i_dom["end_1based"]))
                if start_diff > 10 or end_diff > 10:
                    conflicts.append(
                        {
                            "conflict_code": "PAE_INTERPRO_BOUNDARY_CONFLICT",
                            "source": "domain_caller",
                            "severity": "WARNING",
                            "pae_range": [p_dom["start_1based"], p_dom["end_1based"]],
                            "interpro_range": [i_dom["start_1based"], i_dom["end_1based"]],
                            "interpro_accession": i_dom["accession"],
                            "overlapping_positions": [],
                            "message": (
                                f"PAE domain [{p_dom['start_1based']}, {p_dom['end_1based']}] and "
                                f"InterPro {i_dom['accession']} [{i_dom['start_1based']}, {i_dom['end_1based']}] "
                                f"differ by more than 10 residues (d_start={start_diff}, d_end={end_diff})."
                            ),
                        }
                    )

    # Choose primary domain intervals (prefer curated InterPro domain-type intervals when available, else PAE domains)
    primary_raw = interpro_domains if interpro_domains else pae_domains
    # Merge overlapping domain-type intervals for linker/tail calculation while recording merged provenance
    sorted_intervals = sorted(
        [
            (
                int(d["start_1based"]),
                int(d["end_1based"]),
                str(d["domain_id"]),
                str(d["source"]),
                str(d.get("accession") or d["domain_id"]),
            )
            for d in primary_raw
        ],
        key=lambda x: (x[0], x[1]),
    )
    consensus_domains: list[dict[str, Any]] = []
    for s, e, dom_id, src, acc in sorted_intervals:
        if not consensus_domains or s > consensus_domains[-1]["end_1based"]:
            bfs = [float(ca_bfactors[r]) for r in range(s - 1, e)]
            intra_pae = None
            if pae_matrix is not None:
                p_vals = [
                    float(pae_matrix[r1][r2])
                    for r1 in range(s - 1, e)
                    for r2 in range(s - 1, e)
                ]
                intra_pae = round(sum(p_vals) / len(p_vals), 3) if p_vals else None
            consensus_domains.append(
                {
                    "domain_id": dom_id,
                    "source": src,
                    "start_1based": s,
                    "end_1based": e,
                    "length": e - s + 1,
                    "mean_ca_confidence": round(sum(bfs) / len(bfs), 2) if bfs else 0.0,
                    "min_ca_confidence": round(min(bfs), 2) if bfs else 0.0,
                    "max_ca_confidence": round(max(bfs), 2) if bfs else 0.0,
                    "mean_intra_pae_angstrom": intra_pae,
                    "merged_domain_ids": [dom_id],
                    "merged_accessions": [acc],
                }
            )
        else:
            prev = consensus_domains[-1]
            new_end = max(prev["end_1based"], e)
            prev["end_1based"] = new_end
            prev["length"] = new_end - prev["start_1based"] + 1
            if dom_id not in prev["merged_domain_ids"]:
                prev["merged_domain_ids"].append(dom_id)
            if acc not in prev["merged_accessions"]:
                prev["merged_accessions"].append(acc)
            bfs = [float(ca_bfactors[r]) for r in range(prev["start_1based"] - 1, new_end)]
            prev["mean_ca_confidence"] = round(sum(bfs) / len(bfs), 2) if bfs else 0.0
            prev["min_ca_confidence"] = round(min(bfs), 2) if bfs else 0.0
            prev["max_ca_confidence"] = round(max(bfs), 2) if bfs else 0.0
            if pae_matrix is not None:
                p_vals = [
                    float(pae_matrix[r1][r2])
                    for r1 in range(prev["start_1based"] - 1, new_end)
                    for r2 in range(prev["start_1based"] - 1, new_end)
                ]
                prev["mean_intra_pae_angstrom"] = round(sum(p_vals) / len(p_vals), 3)

    # Compute inter-domain linkers and N/C-terminal tails
    inter_domain_linkers: list[dict[str, Any]] = []
    for idx in range(len(consensus_domains) - 1):
        d1 = consensus_domains[idx]
        d2 = consensus_domains[idx + 1]
        l_start = d1["end_1based"] + 1
        l_end = d2["start_1based"] - 1
        if l_start <= l_end:
            inter_pae = None
            if pae_matrix is not None:
                cross_vals = [
                    (float(pae_matrix[r1][r2]) + float(pae_matrix[r2][r1])) / 2.0
                    for r1 in range(d1["start_1based"] - 1, d1["end_1based"])
                    for r2 in range(d2["start_1based"] - 1, d2["end_1based"])
                ]
                inter_pae = round(sum(cross_vals) / len(cross_vals), 3) if cross_vals else None
            inter_domain_linkers.append(
                {
                    "linker_id": f"LINKER_{idx + 1}_{idx + 2}",
                    "upstream_domain_id": d1["domain_id"],
                    "downstream_domain_id": d2["domain_id"],
                    "start_1based": l_start,
                    "end_1based": l_end,
                    "length": l_end - l_start + 1,
                    "mean_inter_domain_pae_angstrom": inter_pae,
                }
            )

    terminal_regions: list[dict[str, Any]] = []
    if consensus_domains:
        first_start = consensus_domains[0]["start_1based"]
        if first_start > 1:
            terminal_regions.append(
                {
                    "region_type": "n_terminal_tail",
                    "start_1based": 1,
                    "end_1based": first_start - 1,
                    "length": first_start - 1,
                }
            )
        last_end = consensus_domains[-1]["end_1based"]
        if last_end < n_res:
            terminal_regions.append(
                {
                    "region_type": "c_terminal_tail",
                    "start_1based": last_end + 1,
                    "end_1based": n_res,
                    "length": n_res - last_end,
                }
            )

    # Check if recommended_positions fall inside a rigid domain interior core (>3 aa from boundary)
    conflicts.extend(
        detect_domain_core_conflicts(
            consensus_domains=consensus_domains,
            positions=recommended_positions,
        )
    )

    # Update domain_relative_confidence with true per-domain stats
    inter_pae_summary: list[dict[str, Any]] = []
    if pae_matrix is not None and len(consensus_domains) >= 2:
        for i in range(len(consensus_domains)):
            for j in range(i + 1, len(consensus_domains)):
                di = consensus_domains[i]
                dj = consensus_domains[j]
                cross = [
                    (float(pae_matrix[r1][r2]) + float(pae_matrix[r2][r1])) / 2.0
                    for r1 in range(di["start_1based"] - 1, di["end_1based"])
                    for r2 in range(dj["start_1based"] - 1, dj["end_1based"])
                ]
                inter_pae_summary.append(
                    {
                        "domain_a": di["domain_id"],
                        "domain_b": dj["domain_id"],
                        "mean_inter_domain_pae_angstrom": round(sum(cross) / len(cross), 3) if cross else None,
                    }
                )

    domain_rel["status"] = "AVAILABLE" if consensus_domains else "NOT_AVAILABLE"
    domain_rel["domain_segmentation_status"] = "AVAILABLE" if consensus_domains else "NOT_AVAILABLE"
    domain_rel["used_for_candidate_recommendation"] = False
    domain_rel["per_domain_confidence"] = consensus_domains
    domain_rel["inter_domain_pae_summary"] = inter_pae_summary
    domain_rel["reason"] = (
        "Domain-relative confidence and inter-domain PAE blocks are recorded for human review "
        "and are never substituted for local per-residue CA confidence."
    )

    caller_status = "CONFLICT" if conflicts else "AVAILABLE"
    return (
        {
            "status": caller_status,
            "tool": "+".join(sources_used),
            "version": DOMAIN_CALLER_VERSION,
            "sources_used": sources_used,
            "reason": None,
            "pae_parameter_profile": dict(PAE_PARAMETER_PROFILES["v1_2_1_baseline"]),
            "pae_domains": pae_domains,
            "interpro_domains": interpro_domains,
            "interpro_non_domain_annotations": interpro_non_domain_annotations,
            "interpro_annotations": interpro_annotations,
            "consensus_domains": consensus_domains,
            "inter_domain_linkers": inter_domain_linkers,
            "terminal_regions": terminal_regions,
            "conflicts": conflicts,
        },
        domain_rel,
    )


def detect_domain_core_conflicts(
    *,
    consensus_domains: Sequence[Mapping[str, Any]],
    positions: Sequence[int],
) -> list[dict[str, Any]]:
    """Flag positions deeply inside a rigid domain core (>3 residues away from domain boundaries)."""
    pos_set = {int(p) for p in positions}
    conflicts: list[dict[str, Any]] = []
    for dom in consensus_domains:
        start = int(dom["start_1based"])
        end = int(dom["end_1based"])
        length = int(dom["length"])
        mean_conf = float(dom.get("mean_ca_confidence", 0.0))
        if length >= 10 and mean_conf >= 75.0:
            core_hits = [p for p in range(start + 3, end - 2) if p in pos_set]
            if core_hits:
                conflicts.append(
                    {
                        "conflict_code": "CANDIDATE_IN_RIGID_DOMAIN_CORE",
                        "source": "domain_caller",
                        "severity": "WARNING",
                        "domain_id": dom["domain_id"],
                        "domain_range": [start, end],
                        "overlapping_positions": sorted(core_hits),
                        "message": (
                            f"Selected position(s) {sorted(core_hits)} lie inside the interior core "
                            f"of high-confidence domain '{dom['domain_id']}' ({start}-{end}, "
                            f"mean CA confidence {mean_conf:.1f})."
                        ),
                    }
                )
    return conflicts


def analyze_uniprot_tm_signal_annotations(
    *,
    protein_sequence: str,
    protein_sha256: str,
    uniprot_sidecar_data: Mapping[str, Any] | None,
    uniprot_sidecar_sha256: str | None,
    recommended_positions: Sequence[int],
    stage: str = "annotate",
) -> dict[str, Any]:
    """Validate and map local UniProt TM/Signal peptide/Topological/Disulfide/Glycosylation annotations."""
    if uniprot_sidecar_data is None:
        return {
            "status": "NOT_AVAILABLE",
            "tool": None,
            "source": None,
            "version": TM_SIGNAL_CALLER_VERSION,
            "accession": None,
            "entry_version": None,
            "sequence_version": None,
            "sidecar_sha256": None,
            "reason": "No local UniProt annotation sidecar supplied.",
            "features": [],
            "conflicts": [],
        }

    n_res = len(protein_sequence)

    _validate_sidecar_binding(uniprot_sidecar_data, protein_sequence, protein_sha256, label="UniProt", stage=stage)

    # Validate sequence / hash / length against input protein
    if (
        "protein_sha256" in uniprot_sidecar_data
        and uniprot_sidecar_data["protein_sha256"] is not None
        and str(uniprot_sidecar_data["protein_sha256"]) != protein_sha256
    ):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            "UniProt sidecar protein_sha256 does not match input protein sequence.",
            stage=stage,
        )

    seq_block = uniprot_sidecar_data.get("sequence")
    if isinstance(seq_block, Mapping):
        seq_val = seq_block.get("value")
        seq_len = seq_block.get("length")
        if isinstance(seq_val, str):
            clean = "".join(seq_val.upper().split())
            if clean != protein_sequence or sha256_text(clean) != protein_sha256:
                raise NBFoldError(
                    "SIDECAR_MISMATCH",
                    "UniProt sidecar sequence.value does not match input protein sequence.",
                    stage=stage,
                )
        if seq_len is not None and int(seq_len) != n_res:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                f"UniProt sidecar sequence.length ({seq_len}) does not match input protein ({n_res} aa).",
                stage=stage,
            )
    elif isinstance(seq_block, str):
        clean = "".join(seq_block.upper().split())
        if clean != protein_sequence:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "UniProt sidecar sequence string does not match input protein sequence.",
                stage=stage,
            )

    if (
        "protein_length" in uniprot_sidecar_data
        and uniprot_sidecar_data["protein_length"] is not None
        and int(uniprot_sidecar_data["protein_length"]) != n_res
    ):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            f"UniProt sidecar protein_length ({uniprot_sidecar_data['protein_length']}) does not match input protein ({n_res} aa).",
            stage=stage,
        )

    accession = str(
        uniprot_sidecar_data.get("primaryAccession")
        or uniprot_sidecar_data.get("accession")
        or "LOCAL_UNIPROT"
    ).strip()
    audit = uniprot_sidecar_data.get("entryAudit")
    entry_version = (
        int(audit["entryVersion"])
        if isinstance(audit, Mapping) and isinstance(audit.get("entryVersion"), int)
        else uniprot_sidecar_data.get("entry_version")
    )
    seq_version = (
        int(audit["sequenceVersion"])
        if isinstance(audit, Mapping) and isinstance(audit.get("sequenceVersion"), int)
        else uniprot_sidecar_data.get("sequence_version")
    )

    raw_features = uniprot_sidecar_data.get("features")
    if not isinstance(raw_features, list):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            "UniProt sidecar must contain a 'features' list.",
            stage=stage,
        )

    normalized_features: list[dict[str, Any]] = []
    for feat in raw_features:
        if not isinstance(feat, Mapping):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "UniProt feature entry must be a JSON object.",
                stage=stage,
            )
        ftype = str(feat.get("type") or feat.get("feature_type") or "").strip()
        if not ftype:
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                "UniProt feature entry is missing 'type'/'feature_type'.",
                stage=stage,
            )
        desc = str(feat.get("description") or "").strip()

        start: int | None = None
        end: int | None = None
        loc_modifier = "EXACT"
        if isinstance(feat.get("location"), Mapping):
            loc = feat["location"]
            s_obj = loc.get("start")
            e_obj = loc.get("end")
            if isinstance(s_obj, Mapping):
                start = s_obj.get("value")
                if s_obj.get("modifier") and s_obj.get("modifier") != "EXACT":
                    loc_modifier = str(s_obj.get("modifier"))
            elif isinstance(s_obj, int):
                start = s_obj
            if isinstance(e_obj, Mapping):
                end = e_obj.get("value")
                if e_obj.get("modifier") and e_obj.get("modifier") != "EXACT":
                    loc_modifier = str(e_obj.get("modifier"))
            elif isinstance(e_obj, int):
                end = e_obj
        else:
            start = feat.get("start_1based", feat.get("start"))
            end = feat.get("end_1based", feat.get("end"))

        if (
            type(start) is not int
            or type(end) is not int
            or not (1 <= start <= end <= n_res)
        ):
            raise NBFoldError(
                "SIDECAR_MISMATCH",
                f"UniProt feature '{ftype}' coordinates [{start}, {end}] are outside 1..{n_res}.",
                stage=stage,
            )

        evidences: list[str] = []
        raw_evs = feat.get("evidences") or feat.get("evidence_codes") or []
        if isinstance(raw_evs, list):
            for ev in raw_evs:
                if isinstance(ev, Mapping) and ev.get("evidenceCode"):
                    evidences.append(str(ev["evidenceCode"]))
                elif isinstance(ev, str):
                    evidences.append(ev)

        desc_lower = desc.lower()
        if loc_modifier != "EXACT" or any(
            kw in desc_lower for kw in ("probable", "potential", "by similarity", "uncertain")
        ):
            uncertainty_class = "qualifier_or_boundary_uncertain"
        elif "ECO:0000269" in evidences:
            uncertainty_class = "experimental_evidence"
        elif any(c in evidences for c in ("ECO:0000250", "ECO:0000305")):
            uncertainty_class = "curated_or_similarity_evidence"
        elif any(c in evidences for c in ("ECO:0000255", "ECO:0000256", "ECO:0007829")):
            uncertainty_class = "predicted_sequence_model_evidence"
        else:
            uncertainty_class = "unspecified_evidence"

        normalized_features.append(
            {
                "source": "uniprot_sidecar",
                "accession": accession,
                "entry_version": entry_version,
                "sequence_version": seq_version,
                "feature_type": ftype,
                "description": desc,
                "start_1based": start,
                "end_1based": end,
                "length": end - start + 1,
                "evidence_codes": evidences,
                "location_modifier": loc_modifier,
                "uncertainty": uncertainty_class,
            }
        )

    conflicts = detect_uniprot_feature_conflicts(
        features=normalized_features,
        positions=recommended_positions,
    )
    status = "CONFLICT" if conflicts else "AVAILABLE"

    return {
        "status": status,
        "tool": "uniprot_sidecar_annotator",
        "source": "uniprot_sidecar",
        "version": TM_SIGNAL_CALLER_VERSION,
        "accession": accession,
        "entry_version": entry_version,
        "sequence_version": seq_version,
        "sidecar_sha256": uniprot_sidecar_sha256,
        "reason": None,
        "features": normalized_features,
        "conflicts": conflicts,
    }


def detect_uniprot_feature_conflicts(
    *,
    features: Sequence[Mapping[str, Any]],
    positions: Sequence[int],
) -> list[dict[str, Any]]:
    """Detect internal Signal peptide / TM overlaps and candidate position overlaps with protected features."""
    pos_set = {int(p) for p in positions}
    conflicts: list[dict[str, Any]] = []

    sp_ranges: list[tuple[int, int]] = []
    tm_ranges: list[tuple[int, int]] = []
    for f in features:
        ft = str(f["feature_type"]).lower()
        s = int(f["start_1based"])
        e = int(f["end_1based"])
        if ft == "signal peptide" or ft == "signal":
            sp_ranges.append((s, e))
        elif ft in {"transmembrane", "intramembrane"}:
            tm_ranges.append((s, e))

    for sp_s, sp_e in sp_ranges:
        for tm_s, tm_e in tm_ranges:
            overlap = sorted(set(range(sp_s, sp_e + 1)) & set(range(tm_s, tm_e + 1)))
            if overlap:
                conflicts.append(
                    {
                        "conflict_code": "SIGNAL_PEPTIDE_TM_OVERLAP_CONFLICT",
                        "source": "tm_signal_caller",
                        "severity": "ERROR",
                        "signal_peptide_range": [sp_s, sp_e],
                        "transmembrane_range": [tm_s, tm_e],
                        "overlapping_positions": overlap,
                        "message": (
                            f"UniProt sidecar defines overlapping Signal peptide [{sp_s}, {sp_e}] "
                            f"and Transmembrane segment [{tm_s}, {tm_e}] at positions {overlap}."
                        ),
                    }
                )

    for f in features:
        ft = str(f["feature_type"])
        ft_low = ft.lower()
        s = int(f["start_1based"])
        e = int(f["end_1based"])
        if ft_low in {"signal peptide", "signal"}:
            hits = sorted(p for p in range(s, e + 1) if p in pos_set)
            if hits:
                conflicts.append(
                    {
                        "conflict_code": "CANDIDATE_OVERLAPS_SIGNAL_PEPTIDE",
                        "source": "tm_signal_caller",
                        "severity": "WARNING",
                        "feature_type": ft,
                        "feature_range": [s, e],
                        "overlapping_positions": hits,
                        "message": (
                            f"Selected position(s) {hits} overlap annotated Signal peptide [{s}, {e}]."
                        ),
                    }
                )
        elif ft_low in {"transmembrane", "intramembrane"}:
            hits = sorted(p for p in range(s, e + 1) if p in pos_set)
            if hits:
                conflicts.append(
                    {
                        "conflict_code": "CANDIDATE_OVERLAPS_TRANSMEMBRANE_HELIX",
                        "source": "tm_signal_caller",
                        "severity": "WARNING",
                        "feature_type": ft,
                        "feature_range": [s, e],
                        "overlapping_positions": hits,
                        "message": (
                            f"Selected position(s) {hits} overlap annotated {ft} segment [{s}, {e}]."
                        ),
                    }
                )
        elif ft_low in {"disulfide bond", "glycosylation", "cross-link"}:
            # Disulfide bond connects endpoints s and e (or single site s==e for glycosylation)
            critical_sites = {s, e} if ft_low == "disulfide bond" else set(range(s, e + 1))
            hits = sorted(p for p in critical_sites if p in pos_set)
            if hits:
                conflicts.append(
                    {
                        "conflict_code": "CANDIDATE_OVERLAPS_DISULFIDE_OR_GLYCOSYLATION",
                        "source": "tm_signal_caller",
                        "severity": "WARNING",
                        "feature_type": ft,
                        "feature_range": [s, e],
                        "overlapping_positions": hits,
                        "message": (
                            f"Selected position(s) {hits} overlap annotated {ft} site(s) {sorted(critical_sites)}."
                        ),
                    }
                )

    return conflicts


def evaluate_positions_against_structural_annotations(
    *,
    positions: Sequence[int],
    domain_caller_report: Mapping[str, Any],
    tm_signal_caller_report: Mapping[str, Any],
    backbone_ss_report: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Evaluate any set of 1-based positions (recommended or user-reviewed) for structural conflicts."""
    all_conflicts: list[dict[str, Any]] = []

    # 1. Include position-independent cross-source conflicts (e.g. PAE vs InterPro boundary conflict, SP/TM overlap)
    for c in domain_caller_report.get("conflicts", []) or []:
        if c.get("conflict_code") in {
            "PAE_INTERPRO_DOMAIN_COUNT_MISMATCH",
            "PAE_INTERPRO_BOUNDARY_CONFLICT",
        }:
            all_conflicts.append(dict(c))
    for c in tm_signal_caller_report.get("conflicts", []) or []:
        if c.get("conflict_code") == "SIGNAL_PEPTIDE_TM_OVERLAP_CONFLICT":
            all_conflicts.append(dict(c))

    # 2. Evaluate position-dependent conflicts against consensus domains, TM/SP features, and H/E cores
    all_conflicts.extend(
        detect_domain_core_conflicts(
            consensus_domains=domain_caller_report.get("consensus_domains", []) or [],
            positions=positions,
        )
    )
    all_conflicts.extend(
        detect_uniprot_feature_conflicts(
            features=tm_signal_caller_report.get("features", []) or [],
            positions=positions,
        )
    )
    # Deduplicate SIGNAL_PEPTIDE_TM_OVERLAP_CONFLICT if already added
    seen_keys: set[str] = set()
    deduped: list[dict[str, Any]] = []
    for item in all_conflicts:
        key = f"{item.get('conflict_code')}:{item.get('overlapping_positions')}:{item.get('pae_range')}:{item.get('feature_range')}"
        if key not in seen_keys:
            seen_keys.add(key)
            deduped.append(item)

    deduped.extend(
        detect_secondary_structure_conflicts(
            segments=backbone_ss_report.get("segments", []) or [],
            positions=positions,
        )
    )
    return deduped


def build_structural_review_summary(
    *,
    protein_length: int,
    candidate_regions: Sequence[Mapping[str, Any]],
    recommended_positions: Sequence[int],
    domain_caller_report: Mapping[str, Any],
    tm_signal_caller_report: Mapping[str, Any],
    backbone_ss_report: Mapping[str, Any],
) -> dict[str, Any]:
    """Build composite per-residue and per-candidate-region structural context and conflict summary."""
    # Map each 1-based residue to its domain context
    res_domain_ctx: dict[int, str] = {i: "unannotated" for i in range(1, protein_length + 1)}
    for tail in domain_caller_report.get("terminal_regions", []) or []:
        for p in range(int(tail["start_1based"]), int(tail["end_1based"]) + 1):
            res_domain_ctx[p] = str(tail["region_type"])
    for linker in domain_caller_report.get("inter_domain_linkers", []) or []:
        for p in range(int(linker["start_1based"]), int(linker["end_1based"]) + 1):
            res_domain_ctx[p] = f"inter_domain_linker:{linker['linker_id']}"
    for dom in domain_caller_report.get("consensus_domains", []) or []:
        s = int(dom["start_1based"])
        e = int(dom["end_1based"])
        did = str(dom["domain_id"])
        for p in range(s, e + 1):
            if p <= s + 2 or p >= e - 2:
                res_domain_ctx[p] = f"domain_boundary:{did}"
            else:
                res_domain_ctx[p] = f"rigid_domain_core:{did}"

    # Map each 1-based residue to UniProt features
    res_uniprot: dict[int, list[str]] = {i: [] for i in range(1, protein_length + 1)}
    for feat in tm_signal_caller_report.get("features", []) or []:
        ft = str(feat["feature_type"])
        s = int(feat["start_1based"])
        e = int(feat["end_1based"])
        if ft.lower() == "disulfide bond":
            res_uniprot[s].append(f"{ft}({s}-{e})")
            if e != s:
                res_uniprot[e].append(f"{ft}({s}-{e})")
        else:
            for p in range(s, e + 1):
                res_uniprot[p].append(ft)

    # Map each 1-based residue to secondary structure state
    res_ss: dict[int, str] = {i: "NOT_AVAILABLE" for i in range(1, protein_length + 1)}
    for item in backbone_ss_report.get("residue_assignments", []) or []:
        res_ss[int(item["protein_position"])] = str(item["ss_state"])

    candidate_region_contexts: list[dict[str, Any]] = []
    for reg in candidate_regions:
        s = int(reg["start_1based"])
        e = int(reg["end_1based"])
        dom_labels = sorted({res_domain_ctx[p] for p in range(s, e + 1)})
        up_labels = sorted({lbl for p in range(s, e + 1) for lbl in res_uniprot[p]})
        ss_labels = sorted({res_ss[p] for p in range(s, e + 1)})
        candidate_region_contexts.append(
            {
                "start_1based": s,
                "end_1based": e,
                "length": int(reg["length"]),
                "domain_contexts": dom_labels,
                "uniprot_features": up_labels,
                "secondary_structure_states": ss_labels,
            }
        )

    annotation_conflicts = evaluate_positions_against_structural_annotations(
        positions=recommended_positions,
        domain_caller_report=domain_caller_report,
        tm_signal_caller_report=tm_signal_caller_report,
        backbone_ss_report=backbone_ss_report,
    )

    return {
        "schema_version": "nbfold-structural-review-v1",
        "sources_summary": {
            "domain_caller": {
                "status": domain_caller_report["status"],
                "tool": domain_caller_report["tool"],
                "version": domain_caller_report["version"],
                "sources_used": domain_caller_report["sources_used"],
            },
            "tm_signal_caller": {
                "status": tm_signal_caller_report["status"],
                "tool": tm_signal_caller_report["tool"],
                "version": tm_signal_caller_report["version"],
                "accession": tm_signal_caller_report["accession"],
            },
            "backbone_secondary_structure": {
                "status": backbone_ss_report["status"],
                "method": backbone_ss_report["method"],
                "version": backbone_ss_report["version"],
            },
            "external_dssp_secondary_structure": {
                "status": "NOT_AVAILABLE",
                "tool": None,
            },
        },
        "has_conflicts": bool(annotation_conflicts),
        "conflict_count": len(annotation_conflicts),
        "annotation_conflicts": annotation_conflicts,
        "candidate_region_contexts": candidate_region_contexts,
        "requires_explicit_conflict_acknowledgement": bool(annotation_conflicts),
    }


def compare_backbone_ss_with_dssp_reference(
    *,
    backbone_ss_report: Mapping[str, Any],
    dssp_3state_reference: str,
) -> dict[str, Any]:
    """Compare local `backbone_secondary_structure` (`H`/`E`/`C`/`U`) against a 1:1 3-state DSSP reference (`H`/`E`/`C`)."""
    assignments = backbone_ss_report.get("residue_assignments", []) or []
    n_res = len(dssp_3state_reference)
    if len(assignments) != n_res:
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            f"DSSP reference length ({n_res}) does not match backbone_secondary_structure assignment count ({len(assignments)}).",
            stage="benchmark",
        )
    ref_states = list(dssp_3state_reference.strip().upper())
    if any(ch not in {"H", "E", "C"} for ch in ref_states):
        raise NBFoldError(
            "SIDECAR_MISMATCH",
            "DSSP 3-state reference string must contain only 'H', 'E', or 'C'.",
            stage="benchmark",
        )

    pred_states = [str(item["ss_state"]) for item in assignments]
    unassigned_count = sum(st == "U" for st in pred_states)
    assigned_pairs = [
        (r, p) for r, p in zip(ref_states, pred_states) if p in {"H", "E", "C"}
    ]
    assigned_count = len(assigned_pairs)
    exact_matches = sum(1 for r, p in assigned_pairs if r == p)
    q3_accuracy = round(exact_matches / assigned_count, 4) if assigned_count > 0 else 0.0

    confusion: dict[str, dict[str, int]] = {
        r: {"H": 0, "E": 0, "C": 0, "U": 0} for r in ("H", "E", "C")
    }
    for r, p in zip(ref_states, pred_states):
        confusion[r][p if p in {"H", "E", "C", "U"} else "U"] += 1

    per_state: dict[str, dict[str, float | int]] = {}
    for st in ("H", "E", "C"):
        tp = confusion[st][st]
        fp = sum(confusion[other][st] for other in ("H", "E", "C") if other != st)
        fn = sum(confusion[st][other] for other in ("H", "E", "C", "U") if other != st)
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2.0 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0
        per_state[st] = {
            "support": sum(1 for r in ref_states if r == st),
            "predicted_count": sum(1 for p in pred_states if p == st),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": round(prec, 4),
            "recall": round(rec, 4),
            "f1": round(f1, 4),
        }

    def _structured_segments(seq_states: Sequence[str]) -> list[tuple[str, int, int]]:
        segs: list[tuple[str, int, int]] = []
        i = 0
        while i < len(seq_states):
            st = seq_states[i]
            j = i
            while j < len(seq_states) and seq_states[j] == st:
                j += 1
            if st in {"H", "E"}:
                segs.append((st, i + 1, j))
            i = j
        return segs

    ref_segs = _structured_segments(ref_states)
    pred_segs = _structured_segments(pred_states)
    boundary_shifts: list[int] = []
    matched_segments: list[dict[str, Any]] = []
    for r_st, r_s, r_e in ref_segs:
        overlapping = [
            (p_st, p_s, p_e)
            for p_st, p_s, p_e in pred_segs
            if p_st == r_st and max(r_s, p_s) <= min(r_e, p_e)
        ]
        if overlapping:
            best = max(
                overlapping,
                key=lambda x: min(r_e, x[2]) - max(r_s, x[1]) + 1,
            )
            ds = abs(r_s - best[1])
            de = abs(r_e - best[2])
            boundary_shifts.extend([ds, de])
            matched_segments.append(
                {
                    "ss_state": r_st,
                    "dssp_range": [r_s, r_e],
                    "observer_range": [best[1], best[2]],
                    "start_shift_aa": ds,
                    "end_shift_aa": de,
                }
            )

    mean_boundary_shift = (
        round(sum(boundary_shifts) / len(boundary_shifts), 3)
        if boundary_shifts
        else 0.0
    )
    max_boundary_shift = max(boundary_shifts) if boundary_shifts else 0

    mismatches = [
        {
            "position_1based": idx + 1,
            "dssp_state": r,
            "observer_state": p,
            "uncertainty": assignments[idx].get("uncertainty"),
        }
        for idx, (r, p) in enumerate(zip(ref_states, pred_states))
        if r != p
    ]

    return {
        "total_residues": n_res,
        "assigned_residues": assigned_count,
        "unassigned_u_count": unassigned_count,
        "unassigned_u_fraction": round(unassigned_count / n_res, 4) if n_res > 0 else 0.0,
        "q3_accuracy": q3_accuracy,
        "per_state_metrics": per_state,
        "confusion_matrix": confusion,
        "matched_structured_segments": matched_segments,
        "mean_boundary_shift_aa": mean_boundary_shift,
        "max_boundary_shift_aa": max_boundary_shift,
        "mismatch_positions": mismatches,
    }


def compare_pae_parameter_sensitivity(
    *,
    pae_matrix: Sequence[Sequence[float]],
    ca_bfactors: Sequence[float],
    interpro_domains: Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Compare PAE domain segmentation across `v1_2_1_baseline` and `prd_v1_0_0_candidate` profiles."""
    n_res = len(pae_matrix)
    interpro_list = [dict(d) for d in (interpro_domains or []) if d.get("used_for_domain_consensus", True)]
    profile_reports: dict[str, dict[str, Any]] = {}

    for prof_id, prof in PAE_PARAMETER_PROFILES.items():
        doms = _segment_pae_domains(
            pae_matrix,
            ca_bfactors,
            distance_cutoff=float(prof["distance_cutoff_angstrom"]),
            merge_cutoff=float(prof["merge_cutoff_angstrom"]),
            min_subdomain_aa=int(prof["min_subdomain_aa"]),
            min_domain_aa=int(prof["min_domain_aa"]),
            min_protein_aa=int(prof["min_protein_aa"]),
        )
        covered: set[int] = set()
        for d in doms:
            covered.update(range(int(d["start_1based"]), int(d["end_1based"]) + 1))
        unassigned_count = n_res - len(covered)
        linkers: list[dict[str, Any]] = []
        for idx in range(len(doms) - 1):
            l_s = int(doms[idx]["end_1based"]) + 1
            l_e = int(doms[idx + 1]["start_1based"]) - 1
            if l_s <= l_e:
                linkers.append(
                    {
                        "linker_id": f"PAE_LINKER_{idx + 1}_{idx + 2}",
                        "start_1based": l_s,
                        "end_1based": l_e,
                        "length": l_e - l_s + 1,
                    }
                )
        conflicts: list[dict[str, Any]] = []
        if interpro_list:
            if len(doms) != len(interpro_list):
                conflicts.append(
                    {
                        "conflict_code": "PAE_INTERPRO_DOMAIN_COUNT_MISMATCH",
                        "pae_domain_count": len(doms),
                        "interpro_domain_count": len(interpro_list),
                    }
                )
            else:
                for p_dom, i_dom in zip(doms, interpro_list):
                    ds = abs(int(p_dom["start_1based"]) - int(i_dom["start_1based"]))
                    de = abs(int(p_dom["end_1based"]) - int(i_dom["end_1based"]))
                    if ds > 10 or de > 10:
                        conflicts.append(
                            {
                                "conflict_code": "PAE_INTERPRO_BOUNDARY_CONFLICT",
                                "pae_range": [p_dom["start_1based"], p_dom["end_1based"]],
                                "interpro_range": [i_dom["start_1based"], i_dom["end_1based"]],
                                "start_diff_aa": ds,
                                "end_diff_aa": de,
                            }
                        )
        profile_reports[prof_id] = {
            "profile": dict(prof),
            "domain_count": len(doms),
            "domains": doms,
            "linker_intervals": linkers,
            "assigned_residues": len(covered),
            "unassigned_residues": unassigned_count,
            "unassigned_residue_fraction": round(unassigned_count / n_res, 4) if n_res > 0 else 0.0,
            "interpro_conflicts": conflicts,
        }

    base_doms = profile_reports["v1_2_1_baseline"]["domains"]
    cand_doms = profile_reports["prd_v1_0_0_candidate"]["domains"]
    boundary_deltas: list[dict[str, Any]] = []
    if len(base_doms) == len(cand_doms):
        for idx, (b_d, c_d) in enumerate(zip(base_doms, cand_doms), start=1):
            boundary_deltas.append(
                {
                    "domain_index": idx,
                    "baseline_range": [b_d["start_1based"], b_d["end_1based"]],
                    "candidate_range": [c_d["start_1based"], c_d["end_1based"]],
                    "start_shift_aa": abs(int(b_d["start_1based"]) - int(c_d["start_1based"])),
                    "end_shift_aa": abs(int(b_d["end_1based"]) - int(c_d["end_1based"])),
                }
            )

    return {
        "protein_length": n_res,
        "profiles": profile_reports,
        "domain_count_delta": (
            profile_reports["prd_v1_0_0_candidate"]["domain_count"]
            - profile_reports["v1_2_1_baseline"]["domain_count"]
        ),
        "unassigned_fraction_delta": round(
            profile_reports["prd_v1_0_0_candidate"]["unassigned_residue_fraction"]
            - profile_reports["v1_2_1_baseline"]["unassigned_residue_fraction"],
            4,
        ),
        "boundary_deltas": boundary_deltas,
    }

