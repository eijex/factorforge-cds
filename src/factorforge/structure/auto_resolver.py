"""Local uploaded-structure adapter. No remote inference, private cache, or synthetic fallback."""
from __future__ import annotations
from io import StringIO
import hashlib
import math
from typing import Any
import warnings
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1
from factorforge.structure.auto_annotator import recommend_structure_pause_regions

def extract_ca_trace(
    structure_content: str,
    fmt: str,
    chain_id: str = "A",
    model_number: int = 1,
) -> list[dict[str, Any]]:
    """Extract centered 3D C-alpha coordinates and B-factor/pLDDT for the browser 3D viewer."""
    parser = (
        PDBParser(PERMISSIVE=True, QUIET=True)
        if fmt == "pdb"
        else MMCIFParser(QUIET=True)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("trace", StringIO(structure_content))
    models = list(structure.get_models())
    model = models[min(max(0, model_number - 1), len(models) - 1)]
    chain = model[chain_id] if chain_id in model else next(iter(model))
    residues = [r for r in chain.get_residues() if r.id[0] == " " and "CA" in r]
    if not residues:
        return []

    coords = [r["CA"].get_coord() for r in residues]
    cx = float(sum(c[0] for c in coords) / len(coords))
    cy = float(sum(c[1] for c in coords) / len(coords))
    cz = float(sum(c[2] for c in coords) / len(coords))

    trace: list[dict[str, Any]] = []
    for idx, r in enumerate(residues, start=1):
        ca = r["CA"].get_coord()
        bf = float(r["CA"].get_bfactor())
        aa = seq1(r.resname, undef_code="X")
        rx = round(float(ca[0]) - cx, 2)
        ry = round(float(ca[1]) - cy, 2)
        rz = round(float(ca[2]) - cz, 2)
        trace.append(
            {
                "pos": idx,
                "aa": aa,
                "x": rx,
                "y": ry,
                "z": rz,
                "ca": [rx, ry, rz],
                "plddt": round(bf, 1),
            }
        )
    return trace


def build_cartoon_pdb_with_helices(
    structure_content: str,
    fmt: str,
    chain_id: str = "A",
    model_number: int = 1,
) -> str:
    """Build a compact PDB with HELIX records and N/CA/C/O/CB atoms for 3Dmol.js cartoon ribbon rendering."""
    from Bio.PDB.Polypeptide import PPBuilder

    parser = (
        PDBParser(PERMISSIVE=True, QUIET=True)
        if fmt == "pdb"
        else MMCIFParser(QUIET=True)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("cartoon", StringIO(structure_content))
    models = list(structure.get_models())
    model = models[min(max(0, model_number - 1), len(models) - 1)]
    chain = model[chain_id] if chain_id in model else next(iter(model))
    residues = [r for r in chain.get_residues() if r.id[0] == " " and "CA" in r]
    if not residues:
        return ""

    lines: list[str] = ["HEADER    UPLOADED STRUCTURE CARTOON"]
    try:
        ppb = PPBuilder()
        peptides = ppb.build_peptides(chain)
        helix_runs: list[tuple[int, int]] = []
        res_to_idx = {id(r): idx for idx, r in enumerate(residues, start=1)}
        for pp in peptides:
            phi_psi = pp.get_phi_psi_list()
            run_start: int | None = None
            prev_idx: int | None = None
            for r_obj, (phi, psi) in zip(pp, phi_psi):
                idx = res_to_idx.get(id(r_obj))
                if idx is None:
                    continue
                is_h = (
                    phi is not None
                    and psi is not None
                    and -115.0 <= math.degrees(phi) <= -25.0
                    and -80.0 <= math.degrees(psi) <= -5.0
                )
                if is_h:
                    if run_start is None:
                        run_start = idx
                    prev_idx = idx
                else:
                    if run_start is not None and prev_idx is not None and (prev_idx - run_start + 1) >= 4:
                        helix_runs.append((run_start, prev_idx))
                    run_start = None
                    prev_idx = None
            if run_start is not None and prev_idx is not None and (prev_idx - run_start + 1) >= 4:
                helix_runs.append((run_start, prev_idx))

        for h_num, (start_i, end_i) in enumerate(helix_runs, start=1):
            r_start = residues[start_i - 1].resname.strip()[:3]
            r_end = residues[end_i - 1].resname.strip()[:3]
            length = end_i - start_i + 1
            lines.append(
                f"HELIX  {h_num:3d} {h_num:3d} {r_start:>3s} A {start_i:4d}  {r_end:>3s} A {end_i:4d}  1                               {length:5d}"
            )
    except Exception:
        pass

    coords = [r["CA"].get_coord() for r in residues]
    cx = float(sum(c[0] for c in coords) / len(coords))
    cy = float(sum(c[1] for c in coords) / len(coords))
    cz = float(sum(c[2] for c in coords) / len(coords))

    atom_serial = 1
    keep_atoms = ("N", "CA", "C", "O", "CB")
    for idx, r in enumerate(residues, start=1):
        res3 = r.resname.strip()[:3].ljust(3)
        ca_bf = float(r["CA"].get_bfactor())
        for aname in keep_atoms:
            if aname not in r:
                continue
            atom = r[aname]
            coord = atom.get_coord()
            elem = aname[0]
            atom_field = f" {aname:<3s}"
            rx = float(coord[0]) - cx
            ry = float(coord[1]) - cy
            rz = float(coord[2]) - cz
            lines.append(
                f"ATOM  {atom_serial:5d} {atom_field:4s} {res3:3s} A{idx:4d}    "
                f"{rx:8.3f}{ry:8.3f}{rz:8.3f}"
                f"  1.00{ca_bf:6.2f}           {elem:>2s}  "
            )
            atom_serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def auto_resolve_structure_for_protein(
    protein: str, uploaded_content: str | None = None,
    uploaded_format: str | None = None, chain_id: str = "A", model_number: int = 1,
) -> dict[str, Any]:
    """Map an explicitly supplied structure locally; recommendations remain unreviewed."""
    clean = "".join(protein.upper().split()).rstrip("*")
    if not uploaded_content:
        raise ValueError("Structure NOT_AVAILABLE: upload an exact-match PDB/mmCIF file. Remote inference and synthetic fallback are disabled.")
    fmt = "pdb" if uploaded_format == "pdb" else "mmcif"
    rec = recommend_structure_pause_regions(clean, uploaded_content, format=fmt,
        chain_id=chain_id, model_number=model_number, is_alphafold=None)
    return {"content": uploaded_content, "format": fmt, "chain": chain_id,
            "model": model_number, "positions": list(rec.recommended_positions),
            "positions_reviewed": False, "source_label": "User-uploaded structure (provider unverified)",
            "protein_sha256": hashlib.sha256(clean.encode()).hexdigest(),
            "structure_sha256": hashlib.sha256(uploaded_content.encode()).hexdigest(),
            "ca_trace": extract_ca_trace(uploaded_content, fmt, chain_id, model_number),
            "cartoon_pdb": build_cartoon_pdb_with_helices(uploaded_content, fmt, chain_id, model_number)}
