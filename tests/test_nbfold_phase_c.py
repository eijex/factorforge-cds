"""Phase C structural annotation and user review strengthening tests (Job 399).

Covers:
1. PAE + InterPro domain boundary & linker context analyzer (`domain_caller`), including:
   - N-terminal / C-terminal disordered tail vs true inter-domain linker distinction
   - Domain-relative PAE confidence (`intra_domain_mean_pae`, `inter_domain_mean_pae`)
   - Cross-validation between PAE domains and InterPro Pfam/CDD sidecar domains
   - Boundary divergence conflict detection (`PAE_INTERPRO_BOUNDARY_DIVERGENCE`)
   - Mismatched InterPro sidecar rejection (`SIDECAR_MISMATCH`)
2. UniProt TM / signal peptide / topological domain / disulfide / glycosylation annotator (`tm_signal_caller`), including:
   - Verified UniProt sidecar intake (`protein_sha256` / `protein_length` / `sequence`)
   - Signal peptide + transmembrane helix overlap detection (`SIGNAL_PEPTIDE_TM_OVERLAP`)
   - Candidate pause position overlap with Signal peptide, TM helix, Disulfide bond, or Glycosylation site
   - Mismatched UniProt sidecar rejection (`SIDECAR_MISMATCH`)
3. Local backbone geometry (`N, CA, C, O`) secondary-structure observer (`backbone_secondary_structure`), including:
   - Full-backbone (`N, CA, C, O`) NeRF peptide coordinates for right-handed alpha-helix (`H`),
     extended beta-strand (`E`), and flexible coil/linker (`C`)
   - `CA`-only structure fallback to `NOT_AVAILABLE` without inventing `H` or `E`
   - Separation from external DSSP (`secondary_structure_dssp.status == "NOT_AVAILABLE"`)
4. Strengthened user review gate (`REVIEW_CONFLICTS_UNACKNOWLEDGED`, exit code 3) and CLI `--interpro-sidecar` / `--uniprot-sidecar` end-to-end execution.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pytest
from Bio.SeqUtils import seq3
from click.testing import CliRunner

from factorforge.cli.main import cli as cli_main
from factorforge.nbfold.contracts import (
    REVIEW_ACKNOWLEDGEMENT_TYPE,
    REVIEW_SCHEMA_VERSION,
    NBFoldError,
)
from factorforge.nbfold.pipeline import run_nbfold_annotate, run_nbfold_optimize

# 28-aa synthetic test protein (>= 20 aa for PAE domain segmentation)
PHASE_C_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLG"
assert len(PHASE_C_PROTEIN) == 28


def _sha256_str(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _v_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _v_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _v_scale(a: tuple[float, float, float], s: float) -> tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def _v_dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _v_cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _v_norm(a: tuple[float, float, float]) -> tuple[float, float, float]:
    n = math.sqrt(_v_dot(a, a))
    return (a[0] / n, a[1] / n, a[2] / n)


def _nerf_place(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
    bond_len: float,
    bond_angle_deg: float,
    dihedral_deg: float,
) -> tuple[float, float, float]:
    """Place atom D given A, B, C, bond length |CD|, angle BCD, and dihedral ABCD."""
    bc = _v_norm(_v_sub(c, b))
    ab = _v_sub(b, a)
    n = _v_norm(_v_cross(ab, bc))
    m = _v_cross(n, bc)

    theta = math.radians(180.0 - bond_angle_deg)
    phi = math.radians(dihedral_deg)

    dx = bond_len * math.cos(theta)
    dy = bond_len * math.sin(theta) * math.cos(phi)
    dz = bond_len * math.sin(theta) * math.sin(phi)

    d = _v_add(
        c,
        _v_add(
            _v_scale(bc, dx),
            _v_add(_v_scale(m, dy), _v_scale(n, dz)),
        ),
    )
    return d


def _build_full_backbone_pdb(
    sequence: str,
    plddts: list[float],
    phi_psi_per_residue: list[tuple[float, float]],
    chain_id: str = "A",
) -> str:
    """Build a standard PDB with N, CA, C, O backbone coordinates from per-residue (phi, psi)."""
    assert len(sequence) == len(plddts) == len(phi_psi_per_residue)
    # Seed virtual atoms prior to N1
    p_prev_c = (-1.329, 0.0, 0.0)
    p_prev_ca = (-2.400, 1.080, 0.0)
    n_pos = (0.0, 0.0, 0.0)

    lines: list[str] = ["HEADER    SYNTHETIC FULL BACKBONE PDB"]
    serial = 1

    for idx, (aa, bval, (phi_deg, psi_deg)) in enumerate(
        zip(sequence, plddts, phi_psi_per_residue, strict=True), start=1
    ):
        res3 = seq3(aa).upper()
        if idx == 1:
            ca_pos = _nerf_place(p_prev_ca, p_prev_c, n_pos, 1.458, 121.7, 180.0)
        else:
            ca_pos = _nerf_place(p_prev_ca, p_prev_c, n_pos, 1.458, 121.7, 180.0)

        c_pos = _nerf_place(p_prev_c, n_pos, ca_pos, 1.525, 111.2, phi_deg)

        # Next N position (used also to place carbonyl O in trans peptide plane)
        next_n_pos = _nerf_place(n_pos, ca_pos, c_pos, 1.329, 116.2, psi_deg)
        o_pos = _nerf_place(n_pos, ca_pos, c_pos, 1.231, 120.8, psi_deg + 180.0)

        for atom_name, elem, coords in (
            ("N", "N", n_pos),
            ("CA", "C", ca_pos),
            ("C", "C", c_pos),
            ("O", "O", o_pos),
        ):
            lines.append(
                f"ATOM  {serial:5d}  {atom_name:<3s} {res3:>3s} {chain_id:1s}{idx:4d}    "
                f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  1.00{bval:6.2f}           {elem:>2s}"
            )
            serial += 1

        p_prev_ca = ca_pos
        p_prev_c = c_pos
        n_pos = next_n_pos

    lines.append("END")
    return "\n".join(lines) + "\n"


def _build_two_domain_pae_matrix(n: int) -> list[list[float]]:
    """Construct a 2-domain PAE matrix for a 28-aa protein:
    - 1..2: N-terminal disordered tail (high PAE = 18.0)
    - 3..12: Domain 1 (low intra-domain PAE = 2.0)
    - 13..16: True inter-domain linker (high PAE = 16.0)
    - 17..26: Domain 2 (low intra-domain PAE = 2.2, high cross-domain PAE = 14.5)
    - 27..28: C-terminal disordered tail (high PAE = 19.0)
    """
    mat: list[list[float]] = []
    for i in range(1, n + 1):
        row: list[float] = []
        for j in range(1, n + 1):
            if i == j:
                row.append(0.5)
            elif 3 <= i <= 12 and 3 <= j <= 12:
                row.append(2.0)
            elif 17 <= i <= 26 and 17 <= j <= 26:
                row.append(2.2)
            elif (3 <= i <= 12 and 17 <= j <= 26) or (17 <= i <= 26 and 3 <= j <= 12):
                row.append(14.5)
            else:
                row.append(16.5)
        mat.append(row)
    return mat


def test_phase_c_full_backbone_secondary_structure_observation(tmp_path: Path) -> None:
    """Verify full-backbone (`N, CA, C, O`) helix, strand, and coil geometry observation."""
    fasta_path = tmp_path / "target.fasta"
    fasta_path.write_text(f">phase_c_target\n{PHASE_C_PROTEIN}\n", encoding="utf-8")

    # Segment 1..2: coil, 3..12: alpha-helix (-57, -47), 13..16: flexible linker (0, 90),
    # 17..24: beta-strand (-120, +125), 25..28: coil
    phi_psi: list[tuple[float, float]] = []
    plddts: list[float] = []
    for pos in range(1, len(PHASE_C_PROTEIN) + 1):
        if 3 <= pos <= 12:
            phi_psi.append((-57.0, -47.0))
            plddts.append(92.0)
        elif 13 <= pos <= 16:
            phi_psi.append((65.0, 35.0))
            plddts.append(52.0)
        elif 17 <= pos <= 24:
            phi_psi.append((-120.0, 125.0))
            plddts.append(89.0)
        else:
            phi_psi.append((70.0, 20.0))
            plddts.append(62.0)

    pdb_path = tmp_path / "full_backbone.pdb"
    pdb_path.write_text(
        _build_full_backbone_pdb(PHASE_C_PROTEIN, plddts, phi_psi),
        encoding="utf-8",
    )

    out_dir = tmp_path / "ann_out"
    res = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        output_dir=out_dir,
        provider_family="af2",
    )
    ann = res["annotation"]
    ss = ann["independent_structural_tools"]["backbone_secondary_structure"]
    dssp = ann["independent_structural_tools"]["secondary_structure_dssp"]

    assert ss["status"] == "AVAILABLE"
    assert ss["complete_backbone_residues"] == len(PHASE_C_PROTEIN)
    assert ss["total_residues"] == len(PHASE_C_PROTEIN)
    assert dssp["status"] == "NOT_AVAILABLE"

    per_res = {r["protein_position"]: r for r in ss["residue_assignments"]}
    # Helical core residues should be called 'H'
    assert any(per_res[p]["ss_state"] == "H" for p in range(4, 12))
    # Beta-strand core residues should be called 'E'
    assert any(per_res[p]["ss_state"] == "E" for p in range(18, 24))
    # Flexible linker residues should be called 'C'
    assert all(per_res[p]["ss_state"] == "C" for p in range(13, 17))

    # Also verify CA-only structure reports NOT_AVAILABLE for backbone_secondary_structure without inventing H or E
    ca_only_lines = [
        line
        for line in pdb_path.read_text(encoding="utf-8").splitlines()
        if not line.startswith("ATOM") or line[12:16].strip() == "CA"
    ]
    ca_only_pdb = tmp_path / "ca_only.pdb"
    ca_only_pdb.write_text("\n".join(ca_only_lines) + "\n", encoding="utf-8")
    res_ca = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=ca_only_pdb,
        output_dir=tmp_path / "ann_ca_only",
        provider_family="af2",
    )["annotation"]
    ss_ca = res_ca["independent_structural_tools"]["backbone_secondary_structure"]
    assert ss_ca["status"] == "NOT_AVAILABLE"
    assert ss_ca["complete_backbone_residues"] == 0


def test_phase_c_domain_caller_pae_and_interpro_agreement_and_conflict(tmp_path: Path) -> None:
    """Verify PAE + InterPro domain segmentation, N/C-tail vs inter-domain linker classification, and boundary conflict."""
    fasta_path = tmp_path / "target.fasta"
    fasta_path.write_text(f">phase_c_target\n{PHASE_C_PROTEIN}\n", encoding="utf-8")

    phi_psi = [(-57.0, -47.0) if (3 <= p <= 12 or 17 <= p <= 26) else (65.0, 35.0) for p in range(1, 29)]
    plddts = [91.0 if (3 <= p <= 12 or 17 <= p <= 26) else 54.0 for p in range(1, 29)]
    pdb_path = tmp_path / "target.pdb"
    pdb_path.write_text(
        _build_full_backbone_pdb(PHASE_C_PROTEIN, plddts, phi_psi),
        encoding="utf-8",
    )

    pae_path = tmp_path / "pae.json"
    pae_path.write_text(
        json.dumps(
            {
                "pae": _build_two_domain_pae_matrix(len(PHASE_C_PROTEIN)),
                "token_chain_ids": ["A"] * len(PHASE_C_PROTEIN),
                "token_res_ids": list(range(1, len(PHASE_C_PROTEIN) + 1)),
            }
        ),
        encoding="utf-8",
    )

    # Case 1: Agreeing InterPro sidecar (domains [3..12] and [17..26])
    interpro_agree = tmp_path / "interpro_agree.json"
    interpro_agree.write_text(
        json.dumps(
            {
                "protein_sha256": _sha256_str(PHASE_C_PROTEIN),
                "protein_length": len(PHASE_C_PROTEIN),
                "source_database": "InterPro",
                "source_version": "103.0",
                "domains": [
                    {
                        "accession": "PF00001",
                        "name": "N_terminal_domain",
                        "database": "Pfam",
                        "start": 3,
                        "end": 12,
                    },
                    {
                        "accession": "PF00002",
                        "name": "C_terminal_domain",
                        "database": "Pfam",
                        "start": 17,
                        "end": 26,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    res_agree = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        pae_sidecar_path=pae_path,
        interpro_sidecar_path=interpro_agree,
        output_dir=tmp_path / "ann_agree",
        provider_family="af2",
    )["annotation"]

    dom_rel = res_agree["domain_relative_confidence"]
    dom_tool = res_agree["independent_structural_tools"]["domain_boundary_assignment"]

    assert dom_rel["status"] == "AVAILABLE"
    assert dom_rel["domain_segmentation_status"] == "AVAILABLE"
    assert len(dom_rel["per_domain_confidence"]) == 2
    assert dom_rel["per_domain_confidence"][0]["mean_intra_pae_angstrom"] < 3.0
    assert dom_rel["inter_domain_pae_summary"][0]["mean_inter_domain_pae_angstrom"] > 12.0

    assert dom_tool["status"] == "AVAILABLE"
    tail_types = [t["region_type"] for t in dom_tool["terminal_regions"]]
    assert "n_terminal_tail" in tail_types
    assert "c_terminal_tail" in tail_types
    assert len(dom_tool["inter_domain_linkers"]) == 1
    assert dom_tool["inter_domain_linkers"][0]["start_1based"] == 13
    assert dom_tool["inter_domain_linkers"][0]["end_1based"] == 16
    # Recommended positions (13..16) are inside the true inter_domain_linker, so no domain conflict!
    assert res_agree["recommended_positions"] == [13, 14, 15, 16]

    # Case 2: Divergent InterPro sidecar (single wide domain [2..27] spanning the PAE linker)
    interpro_conflict = tmp_path / "interpro_conflict.json"
    interpro_conflict.write_text(
        json.dumps(
            {
                "protein_sha256": _sha256_str(PHASE_C_PROTEIN),
                "protein_length": len(PHASE_C_PROTEIN),
                "source_database": "InterPro",
                "source_version": "103.0",
                "domains": [
                    {
                        "accession": "PF99999",
                        "name": "Merged_Superdomain",
                        "database": "Pfam",
                        "start": 2,
                        "end": 27,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    res_conflict = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        pae_sidecar_path=pae_path,
        interpro_sidecar_path=interpro_conflict,
        output_dir=tmp_path / "ann_conflict",
        provider_family="af2",
    )["annotation"]

    dom_conflict_tool = res_conflict["independent_structural_tools"]["domain_boundary_assignment"]
    assert dom_conflict_tool["status"] == "CONFLICT"
    conflict_codes = [
        c["conflict_code"]
        for c in res_conflict["structural_review_summary"]["annotation_conflicts"]
    ]
    assert "PAE_INTERPRO_DOMAIN_COUNT_MISMATCH" in conflict_codes
    assert "CANDIDATE_IN_RIGID_DOMAIN_CORE" in conflict_codes


def test_phase_c_uniprot_tm_signal_caller_and_review_conflict_gate(tmp_path: Path) -> None:
    """Verify UniProt sidecar annotations, conflict detection, and REVIEW_CONFLICTS_UNACKNOWLEDGED gate."""
    fasta_path = tmp_path / "target.fasta"
    fasta_path.write_text(f">phase_c_target\n{PHASE_C_PROTEIN}\n", encoding="utf-8")

    # Make positions 13..16 low pLDDT so auto_annotator recommends [13, 14, 15, 16]
    phi_psi = [(65.0, 35.0) for _ in range(len(PHASE_C_PROTEIN))]
    plddts = [52.0 if 13 <= p <= 16 else 90.0 for p in range(1, len(PHASE_C_PROTEIN) + 1)]
    pdb_path = tmp_path / "target.pdb"
    pdb_path.write_text(
        _build_full_backbone_pdb(PHASE_C_PROTEIN, plddts, phi_psi),
        encoding="utf-8",
    )

    # UniProt sidecar declaring a Transmembrane helix across 10..22 (overlapping 13..16)
    # and a Signal peptide 1..6 + Disulfide bond at 14
    uniprot_path = tmp_path / "uniprot.json"
    uniprot_path.write_text(
        json.dumps(
            {
                "protein_sha256": _sha256_str(PHASE_C_PROTEIN),
                "protein_length": len(PHASE_C_PROTEIN),
                "accession": "P00001",
                "source_database": "UniProtKB",
                "source_version": "2026_02",
                "features": [
                    {
                        "type": "Signal peptide",
                        "start": 1,
                        "end": 6,
                        "evidences": ["ECO:0000255"],
                    },
                    {
                        "type": "Transmembrane",
                        "start": 10,
                        "end": 22,
                        "description": "Helical",
                        "evidences": ["ECO:0000255"],
                    },
                    {
                        "type": "Disulfide bond",
                        "start": 14,
                        "end": 25,
                        "evidences": ["ECO:0000269"],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    ann_dir = tmp_path / "ann_tm"
    ann_res = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        uniprot_sidecar_path=uniprot_path,
        output_dir=ann_dir,
        provider_family="af2",
    )
    ann = ann_res["annotation"]
    ann_file = Path(ann_res["annotation_path"])

    tm_tool = ann["independent_structural_tools"]["transmembrane_topology"]
    assert tm_tool["status"] == "CONFLICT"
    assert len(tm_tool["features"]) == 3

    rev_summary = ann["structural_review_summary"]
    assert rev_summary["requires_explicit_conflict_acknowledgement"] is True
    codes = {c["conflict_code"] for c in rev_summary["annotation_conflicts"]}
    assert "CANDIDATE_OVERLAPS_TRANSMEMBRANE_HELIX" in codes
    assert "CANDIDATE_OVERLAPS_DISULFIDE_OR_GLYCOSYLATION" in codes

    # Attempt optimize WITHOUT acknowledging conflicts -> must fail with REVIEW_CONFLICTS_UNACKNOWLEDGED (exit code 3)
    rev_unack_path = tmp_path / "review_unack.json"
    base_review = {
        "schema_version": REVIEW_SCHEMA_VERSION,
        "review_type": REVIEW_ACKNOWLEDGEMENT_TYPE,
        "positions_reviewed": True,
        "reviewed_at": "2026-10-13T10:30:00Z",
        "reviewer_label": "phase_c_tester",
        "protein_sha256": ann["protein_sha256"],
        "structure_sha256": ann["structure_sha256"],
        "sidecar_sha256s": ann["sidecar_sha256s"],
        "chain_id": ann["chain_id"],
        "model_number": ann["model_number"],
        "annotation_sha256": ann["annotation_content_sha256"],
        "reviewed_positions": ann["recommended_positions"],
        "config_digest": ann["default_config_digest"],
        "modification_rationale": None,
    }
    rev_unack_path.write_text(json.dumps(base_review), encoding="utf-8")

    with pytest.raises(NBFoldError) as exc_info:
        run_nbfold_optimize(
            protein_fasta_path=fasta_path,
            structure_file_path=pdb_path,
            annotation_path=ann_file,
            review_manifest_path=rev_unack_path,
            uniprot_sidecar_path=uniprot_path,
            output_dir=tmp_path / "opt_unack",
        )
    assert exc_info.value.error_code == "REVIEW_CONFLICTS_UNACKNOWLEDGED"
    assert exc_info.value.exit_code == 3

    # Now provide explicit conflict acknowledgement + rationale -> must succeed!
    rev_ack_path = tmp_path / "review_ack.json"
    ack_review = {
        **base_review,
        "conflicts_acknowledged": True,
        "conflict_resolution_rationale": (
            "Reviewed TM and disulfide overlap at 13..16; intentional pause test for membrane insertion study."
        ),
    }
    rev_ack_path.write_text(json.dumps(ack_review), encoding="utf-8")

    opt_res = run_nbfold_optimize(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        annotation_path=ann_file,
        review_manifest_path=rev_ack_path,
        uniprot_sidecar_path=uniprot_path,
        output_dir=tmp_path / "opt_ack",
    )
    assert opt_res["run_manifest"]["status"] == "COMPLETED"
    assert opt_res["run_manifest"]["review_provenance"]["conflicts_present"] is True
    assert opt_res["run_manifest"]["review_provenance"]["conflicts_acknowledged"] is True
    assert (
        opt_res["observations"]["structural_annotation_summary"]["transmembrane_topology_status"]
        == "CONFLICT"
    )


def test_phase_c_sidecar_mismatch_and_cli_end_to_end(tmp_path: Path) -> None:
    """Verify mismatched InterPro/UniProt sidecar rejection and CLI `--interpro-sidecar` / `--uniprot-sidecar`."""
    fasta_path = tmp_path / "target.fasta"
    fasta_path.write_text(f">phase_c_target\n{PHASE_C_PROTEIN}\n", encoding="utf-8")

    phi_psi = [(65.0, 35.0) for _ in range(len(PHASE_C_PROTEIN))]
    plddts = [52.0 if 13 <= p <= 16 else 90.0 for p in range(1, len(PHASE_C_PROTEIN) + 1)]
    pdb_path = tmp_path / "target.pdb"
    pdb_path.write_text(
        _build_full_backbone_pdb(PHASE_C_PROTEIN, plddts, phi_psi),
        encoding="utf-8",
    )

    bad_uniprot = tmp_path / "bad_uniprot.json"
    bad_uniprot.write_text(
        json.dumps(
            {
                "protein_sha256": "0" * 64,
                "protein_length": len(PHASE_C_PROTEIN),
                "features": [{"type": "Signal peptide", "start": 1, "end": 5}],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(NBFoldError) as exc_info:
        run_nbfold_annotate(
            protein_fasta_path=fasta_path,
            structure_file_path=pdb_path,
            uniprot_sidecar_path=bad_uniprot,
            output_dir=tmp_path / "bad_out",
            provider_family="af2",
        )
    assert exc_info.value.error_code == "SIDECAR_MISMATCH"
    assert exc_info.value.exit_code == 2

    # Valid CLI run with InterPro and UniProt sidecars
    good_interpro = tmp_path / "good_interpro.json"
    good_interpro.write_text(
        json.dumps(
            {
                "protein_sha256": _sha256_str(PHASE_C_PROTEIN),
                "protein_length": len(PHASE_C_PROTEIN),
                "domains": [
                    {"accession": "PF00010", "name": "DomA", "start": 1, "end": 11},
                    {"accession": "PF00011", "name": "DomB", "start": 18, "end": 28},
                ],
            }
        ),
        encoding="utf-8",
    )
    good_uniprot = tmp_path / "good_uniprot.json"
    good_uniprot.write_text(
        json.dumps(
            {
                "protein_sha256": _sha256_str(PHASE_C_PROTEIN),
                "protein_length": len(PHASE_C_PROTEIN),
                "features": [
                    {"type": "Topological domain", "start": 1, "end": 28, "description": "Cytoplasmic"}
                ],
            }
        ),
        encoding="utf-8",
    )

    runner = CliRunner()
    cli_ann_dir = tmp_path / "cli_ann"
    res_cli_ann = runner.invoke(
        cli_main,
        [
            "nbfold",
            "annotate",
            "--protein-fasta",
            str(fasta_path),
            "--structure-file",
            str(pdb_path),
            "--provider-family",
            "af2",
            "--interpro-sidecar",
            str(good_interpro),
            "--uniprot-sidecar",
            str(good_uniprot),
            "--output-dir",
            str(cli_ann_dir),
        ],
    )
    assert res_cli_ann.exit_code == 0, res_cli_ann.output
    stdout_ann = json.loads(res_cli_ann.output.strip())
    assert stdout_ann["conflict_count"] == 0
    assert stdout_ann["requires_explicit_conflict_acknowledgement"] is False


def test_phase_c_sp_tm_overlap_and_ss_core_conflict(tmp_path: Path) -> None:
    """Verify SIGNAL_PEPTIDE_TM_OVERLAP_CONFLICT and CANDIDATE_IN_RIGID_SECONDARY_STRUCTURE_CORE."""
    fasta_path = tmp_path / "target.fasta"
    fasta_path.write_text(f">phase_c_target\n{PHASE_C_PROTEIN}\n", encoding="utf-8")

    # Alpha-helix across 3..12, linker 13..16, beta-strand 17..24
    phi_psi = [
        (-57.0, -47.0)
        if 3 <= p <= 12
        else ((-120.0, 125.0) if 17 <= p <= 24 else (65.0, 35.0))
        for p in range(1, 29)
    ]
    plddts = [52.0 if 13 <= p <= 16 else 92.0 for p in range(1, 29)]
    pdb_path = tmp_path / "target.pdb"
    pdb_path.write_text(
        _build_full_backbone_pdb(PHASE_C_PROTEIN, plddts, phi_psi),
        encoding="utf-8",
    )

    # UniProt sidecar with overlapping Signal peptide [1..8] and Transmembrane [6..20]
    sp_tm_overlap = tmp_path / "sp_tm_overlap.json"
    sp_tm_overlap.write_text(
        json.dumps(
            {
                "protein_sha256": _sha256_str(PHASE_C_PROTEIN),
                "protein_length": len(PHASE_C_PROTEIN),
                "features": [
                    {"type": "Signal peptide", "start": 1, "end": 8},
                    {"type": "Transmembrane", "start": 6, "end": 20},
                ],
            }
        ),
        encoding="utf-8",
    )
    ann_res = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        uniprot_sidecar_path=sp_tm_overlap,
        output_dir=tmp_path / "ann_overlap",
        provider_family="af2",
    )
    ann = ann_res["annotation"]
    codes = {
        c["conflict_code"]
        for c in ann["structural_review_summary"]["annotation_conflicts"]
    }
    assert "SIGNAL_PEPTIDE_TM_OVERLAP_CONFLICT" in codes

    # Now test user manually selecting positions [6, 7, 8] inside the alpha-helix core (3..12)
    # without acknowledging the secondary structure core conflict -> must fail with REVIEW_CONFLICTS_UNACKNOWLEDGED
    ann_no_up_res = run_nbfold_annotate(
        protein_fasta_path=fasta_path,
        structure_file_path=pdb_path,
        output_dir=tmp_path / "ann_no_up",
        provider_family="af2",
    )
    ann2 = ann_no_up_res["annotation"]
    ann2_file = Path(ann_no_up_res["annotation_path"])
    assert ann2["structural_review_summary"]["conflict_count"] == 0

    rev_ss_conflict = tmp_path / "rev_ss_conflict.json"
    rev_ss_conflict.write_text(
        json.dumps(
            {
                "schema_version": REVIEW_SCHEMA_VERSION,
                "review_type": REVIEW_ACKNOWLEDGEMENT_TYPE,
                "positions_reviewed": True,
                "reviewed_at": "2026-10-13T10:35:00Z",
                "reviewer_label": "ss_core_tester",
                "protein_sha256": ann2["protein_sha256"],
                "structure_sha256": ann2["structure_sha256"],
                "sidecar_sha256s": ann2["sidecar_sha256s"],
                "chain_id": ann2["chain_id"],
                "model_number": ann2["model_number"],
                "annotation_sha256": ann2["annotation_content_sha256"],
                "reviewed_positions": [6, 7, 8],
                "config_digest": ann2["default_config_digest"],
                "modification_rationale": "Testing custom positions inside helix core.",
                "conflicts_acknowledged": False,
                "conflict_resolution_rationale": None,
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(NBFoldError) as exc_info:
        run_nbfold_optimize(
            protein_fasta_path=fasta_path,
            structure_file_path=pdb_path,
            annotation_path=ann2_file,
            review_manifest_path=rev_ss_conflict,
            output_dir=tmp_path / "opt_ss_fail",
        )
    assert exc_info.value.error_code == "REVIEW_CONFLICTS_UNACKNOWLEDGED"
    assert (
        "CANDIDATE_IN_RIGID_SECONDARY_STRUCTURE_CORE"
        in exc_info.value.details["conflict_codes"]
    )



@pytest.mark.parametrize("kind", ["interpro", "uniprot"])
@pytest.mark.parametrize("binding", [{}, {"protein_length": "28"}, {"protein_length": 28.5}, {"protein_length": True}, {"sequence": 7}, {"protein_sha256": None}])
def test_phase_c_rejects_missing_or_malformed_binding(kind, binding):
    from factorforge.nbfold.annotations import analyze_domain_annotations, analyze_uniprot_tm_signal_annotations

    kwargs = dict(protein_sequence=PHASE_C_PROTEIN, protein_sha256=_sha256_str(PHASE_C_PROTEIN), recommended_positions=[])
    with pytest.raises(NBFoldError) as exc:
        if kind == "interpro":
            analyze_domain_annotations(**kwargs, ca_bfactors=[90.0] * 28, pae_matrix=None,
                interpro_sidecar_data=dict(binding, entries=[]), interpro_sidecar_sha256="a" * 64,
                base_domain_relative_confidence={})
        else:
            analyze_uniprot_tm_signal_annotations(**kwargs,
                uniprot_sidecar_data=dict(binding, features=[]), uniprot_sidecar_sha256="b" * 64)
    assert exc.value.error_code == "SIDECAR_MISMATCH"


def test_merged_domain_pae_statistics_and_matches_alias():
    from factorforge.nbfold.annotations import analyze_domain_annotations

    matrix = [[float(i + j) for j in range(28)] for i in range(28)]
    report, _ = analyze_domain_annotations(
        protein_sequence=PHASE_C_PROTEIN, protein_sha256=_sha256_str(PHASE_C_PROTEIN),
        ca_bfactors=[90.0] * 28, pae_matrix=matrix, recommended_positions=[],
        interpro_sidecar_sha256="a" * 64, base_domain_relative_confidence={},
        interpro_sidecar_data={"protein_length": 28, "matches": [
            {"accession": "D1", "start": 1, "end": 10},
            {"accession": "D2", "start": 8, "end": 20},
        ]},
    )
    domain = report["consensus_domains"][0]
    assert domain["end_1based"] == 20
    assert domain["mean_intra_pae_angstrom"] == 19.0


def test_degenerate_full_backbone_is_unassigned():
    from types import SimpleNamespace
    from factorforge.nbfold.annotations import analyze_backbone_secondary_structure

    residues = [{name: SimpleNamespace(get_coord=lambda: (0.0, 0.0, 0.0))
                 for name in ("N", "CA", "C", "O")} for _ in range(3)]
    result = analyze_backbone_secondary_structure(polymer_residues=residues,
        protein_sequence="MAA", recommended_positions=[])
    assert result["status"] == "NOT_AVAILABLE"
    assert all(item["ss_state"] == "U" for item in result["residue_assignments"])


@pytest.mark.parametrize("metric", ["score", "evalue"])
def test_interpro_rejects_nonfinite_metrics(metric):
    from factorforge.nbfold.annotations import analyze_domain_annotations

    with pytest.raises(NBFoldError) as exc:
        analyze_domain_annotations(
            protein_sequence=PHASE_C_PROTEIN, protein_sha256=_sha256_str(PHASE_C_PROTEIN),
            ca_bfactors=[90.0] * 28, pae_matrix=None, recommended_positions=[],
            interpro_sidecar_sha256="a" * 64, base_domain_relative_confidence={},
            interpro_sidecar_data={"protein_length": 28, "entries": [
                {"accession": "D1", "start": 1, "end": 10, metric: float("nan")}]},
        )
    assert exc.value.error_code == "SIDECAR_MISMATCH"
