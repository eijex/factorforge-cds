"""Job 400 Phase C specification alignment, InterPro entry-type separation, and synthetic reconstruction benchmark tests.

Covers:
1. PAE parameter profile semantics (`v1_2_1_baseline` vs `prd_v1_0_0_candidate`) and preservation of `v1.2.1` defaults.
2. InterPro `entry_type` separation (`domain` vs `family`, `homologous_superfamily`, `repeat`, `site`) so full-length
   Family/Superfamily spans or short Repeat/Site motifs are preserved in `interpro_non_domain_annotations` without
   collapsing multi-domain boundaries or erasing inter-domain linkers.
3. Public structure benchmark with strictly separated `calibration_set` (`1CRN_A`, `1BPI_A`) and `evaluation_set`
   (`1UBQ_A`, `1PGA_A`, `2DOM_SH3_PDZ`, `2DOM_IG_LC`) comparing:
   - Local `backbone_secondary_structure` (`H`/`E`/`C`/`U`) vs 3-state DSSP references (`H`/`E`/`C`)
   - PAE segmentation sensitivity (`v1_2_1_baseline` vs `prd_v1_0_0_candidate`)
"""

from __future__ import annotations

import hashlib
import io
import json
import math
from typing import Any

import pytest
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq3

from factorforge.nbfold.annotations import (
    PAE_PARAMETER_PROFILES,
    analyze_backbone_secondary_structure,
    analyze_domain_annotations,
    compare_backbone_ss_with_dssp_reference,
    compare_pae_parameter_sensitivity,
)
from factorforge.nbfold.contracts import NBFoldError


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
    bc = _v_norm(_v_sub(c, b))
    ab = _v_sub(b, a)
    n = _v_norm(_v_cross(ab, bc))
    m = _v_cross(n, bc)

    theta = math.radians(180.0 - bond_angle_deg)
    phi = math.radians(dihedral_deg)

    dx = bond_len * math.cos(theta)
    dy = bond_len * math.sin(theta) * math.cos(phi)
    dz = bond_len * math.sin(theta) * math.sin(phi)

    return _v_add(
        c,
        _v_add(
            _v_scale(bc, dx),
            _v_add(_v_scale(m, dy), _v_scale(n, dz)),
        ),
    )


def _build_full_backbone_pdb_from_dssp(
    sequence: str,
    dssp_3state: str,
    plddts: list[float],
    *,
    cap_transitions: dict[int, tuple[float, float]] | None = None,
    chain_id: str = "A",
) -> str:
    """Build a 1:1 full-backbone (N, CA, C, O) PDB trace from sequence and 3-state DSSP with realistic cap transitions."""
    assert len(sequence) == len(dssp_3state) == len(plddts)
    caps = cap_transitions or {}
    p_prev_c = (-1.329, 0.0, 0.0)
    p_prev_ca = (-2.400, 1.080, 0.0)
    n_pos = (0.0, 0.0, 0.0)

    lines: list[str] = ["HEADER    PUBLIC BENCHMARK FULL BACKBONE PDB"]
    serial = 1

    for idx, (aa, st, bval) in enumerate(zip(sequence, dssp_3state, plddts, strict=True), start=1):
        if idx in caps:
            phi_deg, psi_deg = caps[idx]
        elif st == "H":
            phi_deg, psi_deg = (-60.0, -45.0)
        elif st == "E":
            phi_deg, psi_deg = (-122.0, 128.0)
        else:
            phi_deg, psi_deg = (68.0, 24.0)

        res3 = seq3(aa).upper()
        ca_pos = _nerf_place(p_prev_ca, p_prev_c, n_pos, 1.458, 121.7, 180.0)
        c_pos = _nerf_place(p_prev_c, n_pos, ca_pos, 1.525, 111.2, phi_deg)
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


def _parse_polymer_residues(pdb_text: str) -> list[Any]:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("bench", io.StringIO(pdb_text))
    model = next(struct.get_models())
    chain = next(model.get_chains())
    return [res for res in chain.get_residues() if res.id[0] == " "]


def _build_domain_pae_matrix(
    n_res: int,
    domain_spans: list[tuple[int, int]],
    *,
    intra_pae: float = 2.2,
    inter_domain_pae: float = 16.5,
    tail_linker_pae: float = 18.0,
) -> list[list[float]]:
    mat: list[list[float]] = []
    for i in range(1, n_res + 1):
        row: list[float] = []
        i_dom = next((idx for idx, (s, e) in enumerate(domain_spans) if s <= i <= e), None)
        for j in range(1, n_res + 1):
            if i == j:
                row.append(0.5)
                continue
            j_dom = next((idx for idx, (s, e) in enumerate(domain_spans) if s <= j <= e), None)
            if i_dom is not None and i_dom == j_dom:
                row.append(intra_pae)
            elif i_dom is not None and j_dom is not None and i_dom != j_dom:
                row.append(inter_domain_pae)
            else:
                row.append(tail_linker_pae)
        mat.append(row)
    return mat


def build_public_benchmark_report() -> dict[str, Any]:
    """Run separated calibration_set and evaluation_set benchmarks and return a deterministic JSON report."""
    fixtures: list[dict[str, Any]] = [
        # Calibration Set (2 public single-domain structures for diagnostic inspection)
        {
            "dataset_split": "calibration_set",
            "fixture_id": "1CRN_A",
            "pdb_id": "1CRN",
            "chain_id": "A",
            "uniprot_accession": "P01542",
            "protein_name": "Crambin (Crambe hispanica)",
            "sequence": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
            "dssp_3state": "CEEECC" + "H" * 13 + "CCC" + "H" * 8 + "C" + "E" * 4 + "C" * 11,
            # Realistic C-cap dihedral transition at residue 19 (3_10/turn exit)
            "cap_transitions": {19: (-75.0, 45.0)},
            "pae_domain_spans": [(2, 39)],  # 38-aa compact core (< 40 aa threshold in prd_v1_0_0_candidate)
            "inter_domain_pae": 16.5,
            "interpro_entries": [
                {
                    "accession": "IPR001087",
                    "name": "Plant thionin domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 1, "end": 45, "evalue": 1.2e-18, "score": 64.5}],
                },
                {
                    "accession": "IPR036574",
                    "name": "Plant thionin superfamily",
                    "type": "homologous_superfamily",
                    "source_database": "gene3d",
                    "locations": [{"start": 1, "end": 46, "evalue": 3.4e-20, "score": 70.1}],
                },
                {
                    "accession": "IPR018357",
                    "name": "Plant thionin conserved site",
                    "type": "site",
                    "source_database": "prosite",
                    "locations": [{"start": 3, "end": 16, "score": 12.0}],
                },
            ],
        },
        {
            "dataset_split": "calibration_set",
            "fixture_id": "1BPI_A",
            "pdb_id": "1BPI",
            "chain_id": "A",
            "uniprot_accession": "P00974",
            "protein_name": "Bovine pancreatic trypsin inhibitor (BPTI mature chain)",
            "sequence": "RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA",
            "dssp_3state": "CC" + "H" * 4 + "C" * 11 + "E" * 7 + "C" * 4 + "E" * 7 + "C" * 12 + "H" * 8 + "CCC",
            # C-cap on short 3_10 helix (residue 6) makes 3..5 a 3-aa helix; N-cap on strand 18
            "cap_transitions": {6: (-70.0, 40.0), 18: (-65.0, 35.0)},
            "pae_domain_spans": [(3, 55)],  # 53-aa Kunitz domain core
            "inter_domain_pae": 16.5,
            "interpro_entries": [
                {
                    "accession": "IPR002223",
                    "name": "Pancreatic trypsin inhibitor Kunitz domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 4, "end": 55, "evalue": 4.0e-24, "score": 82.0}],
                },
                {
                    "accession": "IPR036880",
                    "name": "Kunitz inhibitor superfamily",
                    "type": "homologous_superfamily",
                    "source_database": "superfamily",
                    "locations": [{"start": 1, "end": 58, "evalue": 1.1e-25, "score": 86.4}],
                },
                {
                    "accession": "IPR020901",
                    "name": "Kunitz inhibitor conserved site",
                    "type": "site",
                    "source_database": "prosite",
                    "locations": [{"start": 30, "end": 51, "score": 15.2}],
                },
            ],
        },
        # Evaluation Set (4 held-out benchmarks: 2 single-domain + 2 two-domain constructs)
        {
            "dataset_split": "evaluation_set",
            "fixture_id": "1UBQ_A",
            "pdb_id": "1UBQ",
            "chain_id": "A",
            "uniprot_accession": "P0CG48",
            "protein_name": "Ubiquitin (Homo sapiens, 1-76)",
            "sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "dssp_3state": (
                "C"
                + "E" * 6
                + "C" * 4
                + "E" * 5
                + "C" * 6
                + "H" * 12
                + "C" * 3
                + "H" * 3
                + "E" * 5
                + "CC"
                + "E" * 2
                + "C" * 7
                + "H" * 3
                + "C" * 6
                + "E" * 6
                + "C" * 5
            ),
            # Realistic C-cap at residue 34 and short 3_10 turn transition at residue 40
            "cap_transitions": {34: (-75.0, 42.0), 40: (-78.0, 48.0)},
            "pae_domain_spans": [(1, 72)],  # 72-aa globular ubiquitin core + 4-aa C-tail (73..76)
            "inter_domain_pae": 16.5,
            "interpro_entries": [
                {
                    "accession": "IPR000626",
                    "name": "Ubiquitin domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 1, "end": 72, "evalue": 2.1e-30, "score": 102.5}],
                },
                {
                    "accession": "IPR029071",
                    "name": "Ubiquitin-like Roll superfamily",
                    "type": "homologous_superfamily",
                    "source_database": "gene3d",
                    "locations": [{"start": 1, "end": 76, "evalue": 1.0e-31, "score": 105.0}],
                },
                {
                    "accession": "IPR019956",
                    "name": "Ubiquitin conserved site",
                    "type": "site",
                    "source_database": "prosite",
                    "locations": [{"start": 26, "end": 44, "score": 18.0}],
                },
            ],
        },
        {
            "dataset_split": "evaluation_set",
            "fixture_id": "1PGA_A",
            "pdb_id": "1PGA",
            "chain_id": "A",
            "uniprot_accession": "P06654",
            "protein_name": "Protein G B1 domain (Streptococcus sp.)",
            "sequence": "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
            "dssp_3state": (
                "C"
                + "E" * 7
                + "C" * 4
                + "E" * 7
                + "C" * 3
                + "H" * 14
                + "C" * 5
                + "E" * 5
                + "C" * 4
                + "E" * 5
                + "C"
            ),
            # Realistic C-cap at residue 36 of the central alpha-helix
            "cap_transitions": {36: (-72.0, 44.0)},
            "pae_domain_spans": [(2, 55)],  # 54-aa B1 domain core
            "inter_domain_pae": 16.5,
            "interpro_entries": [
                {
                    "accession": "IPR015988",
                    "name": "Immunoglobulin G-binding protein G B1 domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 2, "end": 55, "evalue": 5.0e-26, "score": 88.2}],
                },
                {
                    "accession": "IPR014756",
                    "name": "Immunoglobulin E-set superfamily",
                    "type": "homologous_superfamily",
                    "source_database": "superfamily",
                    "locations": [{"start": 1, "end": 56, "evalue": 1.2e-27, "score": 91.0}],
                },
                {
                    "accession": "IPR018097",
                    "name": "IgG-binding B repeat",
                    "type": "repeat",
                    "source_database": "smart",
                    "locations": [{"start": 10, "end": 48, "score": 41.0}],
                },
            ],
        },
        {
            "dataset_split": "evaluation_set",
            "fixture_id": "2DOM_SH3_PDZ",
            "pdb_id": "SYNTH_SH3_PDZ",
            "chain_id": "A",
            "uniprot_accession": "Q9Y6N5_2DOM",
            "protein_name": "Two-domain SH3 (34 aa) + PDZ (47 aa) tandem construct (92 aa)",
            "sequence": (
                "MGS"
                "KELVLALYDYQEKSPREVTMKKGDILTLLNSTNK"
                "GGSGGSG"
                "PRTITVHRKGLGFSIAGGKGSTPYRAGDAGIFVSRIAEGGAAHRAGT"
                "G"
            ),
            "dssp_3state": (
                "CCC"
                + ("E" * 6 + "CC" + "E" * 6 + "CC" + "E" * 6 + "CC" + "E" * 6 + "CCCC")
                + "CCCCCCC"
                + ("E" * 7 + "CCC" + "E" * 6 + "CC" + "H" * 10 + "CCC" + "E" * 7 + "CCC" + "E" * 6)
                + "C"
            ),
            "cap_transitions": {77: (-74.0, 42.0)},
            "pae_domain_spans": [(4, 37), (45, 91)],  # Domain 1 = 34 aa (<40 aa), Linker = 38..44 (7 aa), Domain 2 = 47 aa
            "inter_domain_pae": 16.5,
            "interpro_entries": [
                {
                    "accession": "IPR001452",
                    "name": "SH3 domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 4, "end": 37, "evalue": 1.5e-16, "score": 58.4}],
                },
                {
                    "accession": "IPR001478",
                    "name": "PDZ domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 45, "end": 91, "evalue": 3.2e-21, "score": 74.0}],
                },
                {
                    "accession": "IPR036028",
                    "name": "SH3-PDZ scaffold family",
                    "type": "family",
                    "source_database": "panther",
                    "locations": [{"start": 1, "end": 92, "evalue": 1.0e-40, "score": 140.0}],
                },
            ],
        },
        {
            "dataset_split": "evaluation_set",
            "fixture_id": "2DOM_IG_LC",
            "pdb_id": "3WD5_L_CORE",
            "chain_id": "L",
            "uniprot_accession": "P01834_2DOM",
            "protein_name": "Two-domain Ig light-chain VL (43 aa core) + CL (48 aa core) with semi-rigid elbow PAE (13.5 A)",
            "sequence": (
                "DIQ"
                "MTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKL"
                "RTVAAPSV"
                "FIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESV"
                "ACA"
            ),
            "dssp_3state": (
                "CCC"
                + ("E" * 8 + "CCC" + "E" * 8 + "CCC" + "E" * 9 + "CCC" + "E" * 9)
                + "CCCCCCCC"
                + ("E" * 9 + "CCC" + "E" * 9 + "CCC" + "H" * 6 + "CCC" + "E" * 8 + "CC" + "E" * 5)
                + "CCC"
            ),
            "cap_transitions": {83: (-74.0, 41.0)},
            "pae_domain_spans": [(4, 46), (55, 102)],  # Domain 1 = 43 aa, Linker = 47..54 (8 aa), Domain 2 = 48 aa
            "inter_domain_pae": 13.5,  # Between 12.0 A (v1_2_1_baseline) and 15.0 A (prd_v1_0_0_candidate)
            "interpro_entries": [
                {
                    "accession": "IPR007110",
                    "name": "Immunoglobulin V-set domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 4, "end": 46, "evalue": 2.0e-19, "score": 67.0}],
                },
                {
                    "accession": "IPR003597",
                    "name": "Immunoglobulin C1-set domain",
                    "type": "domain",
                    "source_database": "pfam",
                    "locations": [{"start": 55, "end": 102, "evalue": 8.0e-22, "score": 76.5}],
                },
                {
                    "accession": "IPR036179",
                    "name": "Immunoglobulin-like fold superfamily",
                    "type": "homologous_superfamily",
                    "source_database": "gene3d",
                    "locations": [{"start": 1, "end": 105, "evalue": 1.0e-45, "score": 152.0}],
                },
            ],
        },
    ]

    calibration_results: list[dict[str, Any]] = []
    evaluation_results: list[dict[str, Any]] = []

    for fx in fixtures:
        seq = str(fx["sequence"])
        dssp = str(fx["dssp_3state"])
        assert len(seq) == len(dssp), f"Length mismatch in {fx['fixture_id']}: {len(seq)} vs {len(dssp)}"
        n_res = len(seq)
        spans: list[tuple[int, int]] = fx["pae_domain_spans"]
        plddts = [
            91.0 if any(s <= pos <= e for s, e in spans) else 56.0
            for pos in range(1, n_res + 1)
        ]
        pdb_text = _build_full_backbone_pdb_from_dssp(
            seq,
            dssp,
            plddts,
            cap_transitions=fx.get("cap_transitions"),
            chain_id=str(fx["chain_id"]),
        )
        residues = _parse_polymer_residues(pdb_text)
        ss_report = analyze_backbone_secondary_structure(
            polymer_residues=residues,
            protein_sequence=seq,
            recommended_positions=[],
        )
        dssp_cmp = compare_backbone_ss_with_dssp_reference(
            backbone_ss_report=ss_report,
            dssp_3state_reference=dssp,
        )

        pae_matrix = _build_domain_pae_matrix(
            n_res,
            spans,
            intra_pae=2.2,
            inter_domain_pae=float(fx["inter_domain_pae"]),
            tail_linker_pae=18.0,
        )
        sidecar_payload = {
            "protein_sha256": _sha256_str(seq),
            "protein_length": n_res,
            "entries": fx["interpro_entries"],
        }
        domain_report, _ = analyze_domain_annotations(
            protein_sequence=seq,
            protein_sha256=_sha256_str(seq),
            ca_bfactors=plddts,
            pae_matrix=pae_matrix,
            interpro_sidecar_data=sidecar_payload,
            interpro_sidecar_sha256=_sha256_str(json.dumps(sidecar_payload, sort_keys=True)),
            recommended_positions=[],
            base_domain_relative_confidence={},
        )
        pae_sensitivity = compare_pae_parameter_sensitivity(
            pae_matrix=pae_matrix,
            ca_bfactors=plddts,
            interpro_domains=domain_report["interpro_domains"],
        )

        item_result = {
            "fixture_id": fx["fixture_id"],
            "dataset_split": fx["dataset_split"],
            "pdb_id": fx["pdb_id"],
            "chain_id": fx["chain_id"],
            "uniprot_accession": fx["uniprot_accession"],
            "protein_name": fx["protein_name"],
            "protein_length": n_res,
            "dssp_comparison": dssp_cmp,
            "interpro_entry_separation": {
                "total_annotations": len(domain_report["interpro_annotations"]),
                "domain_entries_used": len(domain_report["interpro_domains"]),
                "non_domain_entries_excluded": len(domain_report["interpro_non_domain_annotations"]),
                "consensus_domains": [
                    [d["start_1based"], d["end_1based"]] for d in domain_report["consensus_domains"]
                ],
                "inter_domain_linkers": [
                    [lk["start_1based"], lk["end_1based"]] for lk in domain_report["inter_domain_linkers"]
                ],
            },
            "pae_parameter_sensitivity": pae_sensitivity,
        }
        if fx["dataset_split"] == "calibration_set":
            calibration_results.append(item_result)
        else:
            evaluation_results.append(item_result)

    def _summarize_split(items: list[dict[str, Any]]) -> dict[str, Any]:
        total_res = sum(int(it["dssp_comparison"]["total_residues"]) for it in items)
        exact_res = sum(
            round(float(it["dssp_comparison"]["q3_accuracy"]) * int(it["dssp_comparison"]["assigned_residues"]))
            for it in items
        )
        mean_q3 = round(
            sum(float(it["dssp_comparison"]["q3_accuracy"]) for it in items) / len(items),
            4,
        )
        overall_q3 = round(exact_res / total_res, 4) if total_res > 0 else 0.0
        max_shift = max(int(it["dssp_comparison"]["max_boundary_shift_aa"]) for it in items)
        base_conflicts = sum(
            len(it["pae_parameter_sensitivity"]["profiles"]["v1_2_1_baseline"]["interpro_conflicts"])
            for it in items
        )
        cand_conflicts = sum(
            len(it["pae_parameter_sensitivity"]["profiles"]["prd_v1_0_0_candidate"]["interpro_conflicts"])
            for it in items
        )
        return {
            "fixture_count": len(items),
            "total_residues": total_res,
            "macro_mean_q3_accuracy": mean_q3,
            "micro_overall_q3_accuracy": overall_q3,
            "max_boundary_shift_aa": max_shift,
            "unassigned_u_total": sum(int(it["dssp_comparison"]["unassigned_u_count"]) for it in items),
            "v1_2_1_baseline_interpro_conflict_fixtures": base_conflicts,
            "prd_v1_0_0_candidate_interpro_conflict_fixtures": cand_conflicts,
        }

    return {
        "schema_version": "nbfold-phase-c-synthetic-benchmark-v1",
        "evidence_class": "SYNTHETIC_RECONSTRUCTION",
        "independent_dssp_executed": False,
        "real_structure_coordinates_used": False,
        "job_id": 400,
        "preserved_default_profile": "v1_2_1_baseline",
        "pae_parameter_profiles": PAE_PARAMETER_PROFILES,
        "calibration_set": {
            "purpose": "Diagnostic parameter & boundary inspection set (never mixed into held-out evaluation metrics).",
            "summary": _summarize_split(calibration_results),
            "fixtures": calibration_results,
        },
        "evaluation_set": {
            "purpose": "Held-out synthetic benchmark set for DSSP Q3 accuracy and PAE sensitivity comparison.",
            "summary": _summarize_split(evaluation_results),
            "fixtures": evaluation_results,
        },
        "recommendation": {
            "preserve_v1_2_1_baseline_as_default": True,
            "rationale": (
                "On constructed synthetic cases in evaluation_set (4 fixtures, 329 residues) and calibration_set (2 fixtures, 104 residues), "
                "v1_2_1_baseline (distance_cutoff=7.0 A, merge_cutoff=12.0 A, min_subdomain_aa=8, min_domain_aa=8) achieves "
                "0 PAE-vs-InterPro domain count conflicts (6/6 fixtures concordant), whereas prd_v1_0_0_candidate "
                "(merge_cutoff=15.0 A, min_subdomain_aa=10, min_domain_aa=40) triggers PAE_INTERPRO_DOMAIN_COUNT_MISMATCH on "
                "3/6 fixtures: it drops compact domains <40 aa (1CRN_A 38 aa, 2DOM_SH3_PDZ SH3 domain 34 aa) and over-merges "
                "two-domain constructs with moderate inter-domain PAE (2DOM_IG_LC 13.5 A), erasing constructed interval gaps."
            ),
        },
    }


def test_pae_parameter_profiles_and_baseline_preservation() -> None:
    """Verify v1.2.1 baseline and PRD v1.0.0 candidate profiles are explicitly defined and v1.2.1 remains default."""
    base = PAE_PARAMETER_PROFILES["v1_2_1_baseline"]
    cand = PAE_PARAMETER_PROFILES["prd_v1_0_0_candidate"]

    assert base["min_protein_aa"] == 20
    assert base["distance_cutoff_angstrom"] == 7.0
    assert base["merge_cutoff_angstrom"] == 12.0
    assert base["min_subdomain_aa"] == 8
    assert base["min_domain_aa"] == 8

    assert cand["min_protein_aa"] == 20
    assert cand["distance_cutoff_angstrom"] == 7.0
    assert cand["merge_cutoff_angstrom"] == 15.0
    assert cand["min_subdomain_aa"] == 10
    assert cand["min_domain_aa"] == 40


def test_interpro_entry_type_separation_preserves_domains_and_linkers() -> None:
    """Verify Family, Homologous_superfamily, Repeat, and Site entries do not collapse multi-domain boundaries."""
    seq = "A" * 92
    sha = _sha256_str(seq)
    pae = _build_domain_pae_matrix(92, [(4, 37), (45, 91)])
    bfactors = [90.0 if (4 <= i <= 37 or 45 <= i <= 91) else 55.0 for i in range(1, 93)]

    sidecar = {
        "protein_sha256": sha,
        "protein_length": 92,
        "entries": [
            {
                "accession": "IPR036028",
                "name": "Full-length scaffold family",
                "type": "family",
                "source_database": "panther",
                "locations": [{"start": 1, "end": 92, "score": 140.0}],
            },
            {
                "accession": "IPR099999",
                "name": "Full-length superfamily",
                "type": "Homologous_superfamily",
                "source_database": "gene3d",
                "locations": [{"start": 1, "end": 92, "score": 135.0}],
            },
            {
                "accession": "IPR001452",
                "name": "SH3 domain",
                "type": "Domain",
                "source_database": "pfam",
                "locations": [{"start": 4, "end": 37, "evalue": 1.2e-15}],
            },
            {
                "accession": "IPR001478",
                "name": "PDZ domain",
                "type": "domain",
                "source_database": "pfam",
                "locations": [{"start": 45, "end": 91, "evalue": 2.4e-20}],
            },
            {
                "accession": "IPR018097",
                "name": "Internal repeat motif",
                "type": "Repeat",
                "source_database": "smart",
                "locations": [{"start": 10, "end": 28, "score": 25.0}],
            },
            {
                "accession": "IPR019956",
                "name": "Conserved binding site",
                "type": "Active_site",
                "source_database": "prosite",
                "locations": [{"start": 50, "end": 58, "score": 14.0}],
            },
        ],
    }

    report, dom_rel = analyze_domain_annotations(
        protein_sequence=seq,
        protein_sha256=sha,
        ca_bfactors=bfactors,
        pae_matrix=pae,
        interpro_sidecar_data=sidecar,
        interpro_sidecar_sha256=_sha256_str(json.dumps(sidecar, sort_keys=True)),
        recommended_positions=[],
        base_domain_relative_confidence={},
    )

    assert report["status"] == "AVAILABLE"
    assert report["pae_parameter_profile"]["profile_id"] == "v1_2_1_baseline"
    assert len(report["interpro_annotations"]) == 6
    assert len(report["interpro_domains"]) == 2
    assert len(report["interpro_non_domain_annotations"]) == 4

    non_domain_types = {item["entry_type"] for item in report["interpro_non_domain_annotations"]}
    assert non_domain_types == {"family", "homologous_superfamily", "repeat", "site"}
    assert all(not item["used_for_domain_consensus"] for item in report["interpro_non_domain_annotations"])
    assert all(item["used_for_domain_consensus"] for item in report["interpro_domains"])

    # Crucially, consensus_domains has 2 distinct domains [4, 37] and [45, 91], NOT a collapsed [1, 92] span
    assert [(d["start_1based"], d["end_1based"]) for d in report["consensus_domains"]] == [
        (4, 37),
        (45, 91),
    ]
    assert len(report["inter_domain_linkers"]) == 1
    linker = report["inter_domain_linkers"][0]
    assert (linker["start_1based"], linker["end_1based"], linker["length"]) == (38, 44, 7)
    assert report["conflicts"] == []
    assert dom_rel["status"] == "AVAILABLE"


def test_interpro_overlapping_domain_entries_record_merged_provenance_and_invalid_type_fails() -> None:
    """Verify overlapping Domain entries record merged provenance and invalid entry_type raises SIDECAR_MISMATCH."""
    seq = "A" * 40
    sha = _sha256_str(seq)
    bfactors = [88.0] * 40

    sidecar = {
        "protein_sha256": sha,
        "protein_length": 40,
        "entries": [
            {
                "accession": "PF00001",
                "type": "domain",
                "source_database": "pfam",
                "locations": [{"start": 4, "end": 22}],
            },
            {
                "accession": "CD00001",
                "type": "domain",
                "source_database": "cdd",
                "locations": [{"start": 6, "end": 25}],
            },
        ],
    }
    report, _ = analyze_domain_annotations(
        protein_sequence=seq,
        protein_sha256=sha,
        ca_bfactors=bfactors,
        pae_matrix=None,
        interpro_sidecar_data=sidecar,
        interpro_sidecar_sha256="sha",
        recommended_positions=[],
        base_domain_relative_confidence={},
    )
    assert len(report["consensus_domains"]) == 1
    cdom = report["consensus_domains"][0]
    assert (cdom["start_1based"], cdom["end_1based"]) == (4, 25)
    assert cdom["merged_domain_ids"] == ["PF00001_1", "CD00001_1"]
    assert cdom["merged_accessions"] == ["PF00001", "CD00001"]

    bad_sidecar = {
        "protein_sha256": sha,
        "protein_length": 40,
        "entries": [
            {
                "accession": "IPR_BAD",
                "type": "unknown_annotation_type",
                "locations": [{"start": 1, "end": 20}],
            }
        ],
    }
    with pytest.raises(NBFoldError) as exc_info:
        analyze_domain_annotations(
            protein_sequence=seq,
            protein_sha256=sha,
            ca_bfactors=bfactors,
            pae_matrix=None,
            interpro_sidecar_data=bad_sidecar,
            interpro_sidecar_sha256="sha",
            recommended_positions=[],
            base_domain_relative_confidence={},
        )
    assert exc_info.value.error_code == "SIDECAR_MISMATCH"


def test_public_benchmark_calibration_and_evaluation_sets_separated() -> None:
    """Verify separated calibration_set and evaluation_set DSSP Q3 accuracy and PAE parameter sensitivity."""
    report = build_public_benchmark_report()

    cal_ids = {f["fixture_id"] for f in report["calibration_set"]["fixtures"]}
    eval_ids = {f["fixture_id"] for f in report["evaluation_set"]["fixtures"]}
    assert cal_ids == {"1CRN_A", "1BPI_A"}
    assert eval_ids == {"1UBQ_A", "1PGA_A", "2DOM_SH3_PDZ", "2DOM_IG_LC"}
    assert cal_ids.isdisjoint(eval_ids)

    cal_sum = report["calibration_set"]["summary"]
    eval_sum = report["evaluation_set"]["summary"]

    # Local backbone secondary structure observer achieves high Q3 agreement (>= 0.90) with 0 unassigned U on full-backbone PDBs
    assert cal_sum["unassigned_u_total"] == 0
    assert eval_sum["unassigned_u_total"] == 0
    assert cal_sum["macro_mean_q3_accuracy"] >= 0.92
    assert eval_sum["macro_mean_q3_accuracy"] >= 0.92
    assert cal_sum["max_boundary_shift_aa"] <= 2
    assert eval_sum["max_boundary_shift_aa"] <= 2

    # PAE sensitivity comparison: v1_2_1_baseline has 0 InterPro domain count conflicts across both sets,
    # whereas prd_v1_0_0_candidate has 1 conflict in calibration_set (1CRN_A < 40 aa) and 2 conflicts in evaluation_set
    # (2DOM_SH3_PDZ 34-aa domain dropped, 2DOM_IG_LC over-merged at merge_cutoff=15.0 A)
    assert cal_sum["v1_2_1_baseline_interpro_conflict_fixtures"] == 0
    assert cal_sum["prd_v1_0_0_candidate_interpro_conflict_fixtures"] == 1
    assert eval_sum["v1_2_1_baseline_interpro_conflict_fixtures"] == 0
    assert eval_sum["prd_v1_0_0_candidate_interpro_conflict_fixtures"] == 2
    assert report["recommendation"]["preserve_v1_2_1_baseline_as_default"] is True
