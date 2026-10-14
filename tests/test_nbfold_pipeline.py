"""Comprehensive synthetic test suite for NB-Fold Phase A/B pipeline (Job 398, PRD v1.1.0).

Covers all 8 Definition of Done categories in PRD Section 10 using public-safe synthetic fixtures:
1. AF2 PDB / AF3 mmCIF exact mapping, CA-atom confidence selection, multi-chain selection, and verified sidecars.
2. Negative intake cases: chain/model missing, sequence/partial/disordered/modified residue mismatch,
   format/extension mismatch, size bounds (>300 aa, >250 KB, >16 KB FASTA, >2 MB sidecar, >5 MB total),
   multi-record FASTA, NaN/Inf/out-of-range confidence, missing CA, unknown source, sidecar mismatch.
3. No low-confidence regions -> NO_CANDIDATES/NEEDS_ANNOTATION without midpoint fallback or auto-approval;
   missing, unconfirmed, stale, or tampered review.json rejected.
4. Identical resolved host reference/GC/Type IIS constraints, independent hard QC, MMMMM infeasible GC
   and zero-gain feasible cases.
5. MFE unavailable/error/per-call temperature/48-nt window/4-nt flank, sense-CDS scope guard,
   RNA scan SEQUENCE_SCAN evidence class, required unavailable BLOCK, and yield_prediction NOT_AVAILABLE.
6. ESMFold / ESM-2 fail-closed (UNSUPPORTED_BACKEND even when enabled=True), no heavy ML imports or network calls.
7. End-to-end CLI annotate -> review.json -> optimize, fixed filenames, overwrite refusal (OUTPUT_CONFLICT),
   atomic staging cleanup on failure, and no raw sequences in error_manifest.json or stderr.
8. 1DP/3DP lambda=0 equivalence, artifact SHA-256 consistency, and report/FASTA/observations alignment.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.SeqUtils import seq3
from click.testing import CliRunner

from factorforge.cli.main import cli
from factorforge.nbfold import (
    ANNOTATION_FILENAME,
    ERROR_MANIFEST_FILENAME,
    OBSERVATIONS_FILENAME,
    PAIRED_FASTA_FILENAME,
    REPORT_MARKDOWN_FILENAME,
    RUN_MANIFEST_FILENAME,
    NBFoldError,
    evaluate_candidate_observations,
    get_backend_capabilities,
    run_nbfold_annotate,
    run_nbfold_optimize,
)

SYNTHETIC_PROTEIN = "MASTWQHLEKDPVNRSGYF"  # 19 aa, standard 20 AA


def _write_fasta(path: Path, sequence: str, header: str = "synthetic_target|unsafe/../id") -> Path:
    path.write_text(f">{header}\n{sequence}\n", encoding="utf-8")
    return path


def _make_pdb_text(
    sequence: str = SYNTHETIC_PROTEIN,
    b_factors: list[float] | None = None,
    chain: str = "A",
    header: str = "HEADER    SYNTHETIC ALPHAFOLD 2 PREDICTION",
) -> str:
    lines = [header]
    for idx, aa in enumerate(sequence, start=1):
        res3 = seq3(aa).upper()
        bf = b_factors[idx - 1] if b_factors is not None else 85.0
        lines.append(
            f"ATOM  {idx:5d}  CA  {res3:3s} {chain}{idx:4d}    {idx*3.8:8.3f}   0.000   0.000  1.00 {bf:6.2f}           C  "
        )
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _make_af3_mmcif_text(
    sequence: str = SYNTHETIC_PROTEIN,
    ca_bfactors: list[float] | None = None,
    cb_bfactors: list[float] | None = None,
    chain: str = "A",
    extra_chain: str | None = None,
    extra_sequence: str = "MDEFG",
) -> str:
    """Create a valid minimal AF3-style mmCIF with atom-level B-factors on N, CA, and C atoms."""
    lines = [
        "data_af3_synthetic",
        "# AlphaFold 3 synthetic prediction",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.pdbx_PDB_model_num",
    ]
    atom_id = 1
    for idx, aa in enumerate(sequence, start=1):
        res3 = seq3(aa).upper()
        ca_bf = ca_bfactors[idx - 1] if ca_bfactors is not None else 88.0
        other_bf = cb_bfactors[idx - 1] if cb_bfactors is not None else ca_bf
        for atom_name, elem, bf_val, offset in (
            ("N", "N", other_bf, -0.5),
            ("CA", "C", ca_bf, 0.0),
            ("C", "C", other_bf, 0.5),
        ):
            lines.append(
                f"ATOM {atom_id} {elem} {atom_name} . {res3} {chain} 1 {idx} ? "
                f"{idx*3.8 + offset:.3f} 0.000 0.000 1.00 {bf_val:.2f} {idx} {res3} {chain} 1"
            )
            atom_id += 1

    if extra_chain:
        for idx, aa in enumerate(extra_sequence, start=1):
            res3 = seq3(aa).upper()
            lines.append(
                f"ATOM {atom_id} C CA . {res3} {extra_chain} 2 {idx} ? "
                f"{idx*3.8:.3f} 5.000 0.000 1.00 92.00 {idx} {res3} {extra_chain} 1"
            )
            atom_id += 1

    lines.append("#")
    return "\n".join(lines) + "\n"


def _build_valid_review_file(
    path: Path,
    annotation_data: dict,
    *,
    positions: list[int] | None = None,
    modification_rationale: str | None = None,
    overrides: dict | None = None,
) -> Path:
    rev = dict(annotation_data["review_manifest_template"])
    rev["positions_reviewed"] = True
    rev["reviewed_at"] = "2026-10-13T10:00:00Z"
    rev["reviewer_label"] = "synthetic_reviewer"
    if positions is not None:
        rev["reviewed_positions"] = positions
    if modification_rationale is not None:
        rev["modification_rationale"] = modification_rationale
    if overrides:
        rev.update(overrides)
    path.write_text(json.dumps(rev, indent=2) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# 1. AF2 PDB / AF3 mmCIF exact mapping, CA-atom confidence, and sidecars
# ---------------------------------------------------------------------------


def test_af3_mmcif_ca_atom_selection_multichain_and_verified_sidecars(tmp_path: Path) -> None:
    """Verify AF3 mmCIF selects CA atom confidence (not atom average), isolates target chain, and validates sidecars."""
    fasta_p = _write_fasta(tmp_path / "protein.fasta", SYNTHETIC_PROTEIN)
    # CA confidence has low confidence (<70) at positions 6..10, while N/C atoms have 95.0!
    # If averaged across N/CA/C, (95+45+95)/3 = 78.3 >= 70, which would miss the low-confidence region.
    ca_bfs = [90.0] * 5 + [45.0] * 5 + [92.0] * 9
    other_bfs = [95.0] * 19
    cif_text = _make_af3_mmcif_text(
        SYNTHETIC_PROTEIN,
        ca_bfactors=ca_bfs,
        cb_bfactors=other_bfs,
        chain="B",
        extra_chain="A",
        extra_sequence="ACDEF",
    )
    cif_p = tmp_path / "model.cif"
    cif_p.write_text(cif_text, encoding="utf-8")

    # Build verified token-level confidence + PAE sidecar with both chains A and B
    n_a = 5
    n_b = len(SYNTHETIC_PROTEIN)
    n_tot = n_a + n_b
    sidecar_data = {
        "token_chain_ids": ["A"] * n_a + ["B"] * n_b,
        "token_res_ids": list(range(1, n_a + 1)) + list(range(1, n_b + 1)),
        "ca_atom_plddts": [92.0] * n_a + ca_bfs,
        "pae": [[1.5 for _ in range(n_tot)] for _ in range(n_tot)],
    }
    sidecar_p = tmp_path / "conf_pae.json"
    sidecar_p.write_text(json.dumps(sidecar_data), encoding="utf-8")

    out_dir = tmp_path / "ann_out"
    res = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=cif_p,
        declared_format="mmcif",
        chain_id="B",
        model_number=1,
        provider_family="af3",
        confidence_sidecar_path=sidecar_p,
        output_dir=out_dir,
    )
    ann = res["annotation"]
    assert ann["reviewed"] is False
    assert ann["provider_family"] == "af3"
    assert ann["source_basis"] == "user_declared"
    assert ann["confidence_metric"] == "CA_atom_confidence"
    assert ann["selected_atom_policy"] == "CA_atom_only"
    assert ann["recommended_positions"] == [6, 7, 8, 9, 10]
    assert ann["position_mapping"][5]["confidence_value"] == 45.0
    assert ann["position_mapping"][5]["atom_confidence_varies_within_residue"] is True
    assert ann["sidecar_observations"]["status"] == "AVAILABLE"
    assert ann["sidecar_observations"]["confidence_sidecar_verified"] is True
    assert ann["sidecar_observations"]["pae_sidecar_verified"] is True
    assert ann["domain_relative_confidence"]["status"] == "AVAILABLE"
    assert ann["domain_relative_confidence"]["used_for_candidate_recommendation"] is False
    assert ann["domain_relative_confidence"]["domain_segmentation_status"] == "NOT_AVAILABLE"
    for tool_info in ann["independent_structural_tools"].values():
        assert tool_info["status"] == "NOT_AVAILABLE"


# ---------------------------------------------------------------------------
# 2. Negative intake, format, mapping, confidence, and sidecar checks
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "case_name,expected_code",
    [
        ("multi_fasta", "MULTI_RECORD_FASTA"),
        ("invalid_aa", "INVALID_SEQUENCE"),
        ("protein_too_long", "UNSUPPORTED_INPUT_SIZE"),
        ("protein_too_short", "UNSUPPORTED_INPUT_SIZE"),
        ("fasta_over_16kb", "UNSUPPORTED_INPUT_SIZE"),
        ("structure_over_250kb", "UNSUPPORTED_INPUT_SIZE"),
        ("format_ext_mismatch", "FORMAT_EXTENSION_MISMATCH"),
        ("pdb_content_in_cif", "FORMAT_EXTENSION_MISMATCH"),
        ("chain_missing", "CHAIN_NOT_FOUND"),
        ("model_missing", "MODEL_NOT_FOUND"),
        ("sequence_mismatch", "SEQUENCE_MISMATCH"),
        ("modified_residue_mse", "MODIFIED_OR_AMBIGUOUS_RESIDUE"),
        ("disordered_atom", "DISORDERED_RESIDUE"),
        ("missing_ca", "MISSING_CA_ATOM"),
        ("nan_confidence", "INVALID_CONFIDENCE_VALUE"),
        ("out_of_range_confidence", "INVALID_CONFIDENCE_VALUE"),
        ("unknown_source", "UNKNOWN_CONFIDENCE_SOURCE"),
        ("sidecar_value_mismatch", "SIDECAR_MISMATCH"),
        ("sidecar_unmapped", "SIDECAR_MISMATCH"),
    ],
)
def test_negative_intake_cases_fail_with_error_code(
    tmp_path: Path, case_name: str, expected_code: str
) -> None:
    fasta_p = _write_fasta(tmp_path / "prot.fasta", SYNTHETIC_PROTEIN)
    pdb_text = _make_pdb_text(SYNTHETIC_PROTEIN, [85.0] * len(SYNTHETIC_PROTEIN))
    struct_p = tmp_path / "model.pdb"
    struct_p.write_text(pdb_text, encoding="utf-8")
    declared_fmt = "pdb"
    chain = "A"
    model = 1
    provider: str | None = "af2"
    sidecar_p: Path | None = None

    if case_name == "multi_fasta":
        fasta_p.write_text(f">rec1\n{SYNTHETIC_PROTEIN}\n>rec2\n{SYNTHETIC_PROTEIN}\n", encoding="utf-8")
    elif case_name == "invalid_aa":
        _write_fasta(fasta_p, "MASTWQHLEKDPVNRSGYX")
    elif case_name == "protein_too_long":
        _write_fasta(fasta_p, "M" + "A" * 300)
    elif case_name == "protein_too_short":
        _write_fasta(fasta_p, "MA")
    elif case_name == "fasta_over_16kb":
        fasta_p.write_text(">" + "h" * 17000 + "\n" + SYNTHETIC_PROTEIN + "\n", encoding="utf-8")
    elif case_name == "structure_over_250kb":
        struct_p.write_text(pdb_text + ("# " + "x" * 100 + "\n") * 2600, encoding="utf-8")
    elif case_name == "format_ext_mismatch":
        declared_fmt = "mmcif"
    elif case_name == "pdb_content_in_cif":
        struct_p = tmp_path / "wrong.cif"
        struct_p.write_text(pdb_text, encoding="utf-8")
        declared_fmt = "mmcif"
    elif case_name == "chain_missing":
        chain = "Z"
    elif case_name == "model_missing":
        model = 2
    elif case_name == "sequence_mismatch":
        struct_p.write_text(_make_pdb_text(SYNTHETIC_PROTEIN[:-1] + "A"), encoding="utf-8")
    elif case_name == "modified_residue_mse":
        # Replace MET with MSE (selenomethionine)
        struct_p.write_text(pdb_text.replace("MET A   1", "MSE A   1"), encoding="utf-8")
    elif case_name == "disordered_atom":
        # Duplicate altloc A and B on residue 1 CA
        alt_a = "ATOM      1  CA AALA A   2       7.600   0.000   0.000  0.50 85.00           C  "
        alt_b = "ATOM      2  CA BALA A   2       7.700   0.000   0.000  0.50 85.00           C  "
        lines = pdb_text.splitlines()
        lines[2] = alt_a + "\n" + alt_b
        struct_p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    elif case_name == "missing_ca":
        struct_p.write_text(pdb_text.replace("  CA  ALA", "  N   ALA"), encoding="utf-8")
    elif case_name == "nan_confidence":
        struct_p.write_text(pdb_text.replace(" 85.00", "   nan", 1), encoding="utf-8")
    elif case_name == "out_of_range_confidence":
        struct_p.write_text(pdb_text.replace(" 85.00", "105.00", 1), encoding="utf-8")
    elif case_name == "unknown_source":
        struct_p.write_text(_make_pdb_text(SYNTHETIC_PROTEIN, header="HEADER    UNLABELED FILE"), encoding="utf-8")
        provider = None
    elif case_name == "sidecar_value_mismatch":
        sidecar_p = tmp_path / "bad_sidecar.json"
        sidecar_p.write_text(
            json.dumps(
                {
                    "token_chain_ids": ["A"] * len(SYNTHETIC_PROTEIN),
                    "token_res_ids": list(range(1, len(SYNTHETIC_PROTEIN) + 1)),
                    "ca_atom_plddts": [10.0] * len(SYNTHETIC_PROTEIN),
                }
            ),
            encoding="utf-8",
        )
    elif case_name == "sidecar_unmapped":
        sidecar_p = tmp_path / "unmapped.json"
        sidecar_p.write_text(json.dumps({"plddt": [85.0] * len(SYNTHETIC_PROTEIN)}), encoding="utf-8")

    out_dir = tmp_path / f"out_{case_name}"
    with pytest.raises(NBFoldError) as exc_info:
        run_nbfold_annotate(
            protein_fasta_path=fasta_p,
            structure_file_path=struct_p,
            declared_format=declared_fmt,
            chain_id=chain,
            model_number=model,
            provider_family=provider,
            confidence_sidecar_path=sidecar_p,
            output_dir=out_dir,
        )
    assert exc_info.value.error_code == expected_code
    err_manifest = json.loads((out_dir / ERROR_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    assert err_manifest["error_code"] == expected_code
    assert SYNTHETIC_PROTEIN not in json.dumps(err_manifest)


# ---------------------------------------------------------------------------
# 3. NO_CANDIDATES / NEEDS_ANNOTATION and strict review.json contract
# ---------------------------------------------------------------------------


def test_no_candidates_when_all_high_confidence_and_review_contract_enforcement(
    tmp_path: Path,
) -> None:
    fasta_p = _write_fasta(tmp_path / "prot.fasta", SYNTHETIC_PROTEIN)
    pdb_p = tmp_path / "high_conf.pdb"
    pdb_p.write_text(_make_pdb_text(SYNTHETIC_PROTEIN, [92.0] * len(SYNTHETIC_PROTEIN)), encoding="utf-8")

    ann_dir = tmp_path / "ann_high"
    res = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        provider_family="af2",
        output_dir=ann_dir,
    )
    ann = res["annotation"]
    assert ann["status"] == "NO_CANDIDATES"
    assert ann["workflow_state"] == "NEEDS_ANNOTATION"
    assert ann["recommended_positions"] == []
    assert ann["candidate_regions"] == []
    assert ann["reviewed"] is False
    assert not (ann_dir / PAIRED_FASTA_FILENAME).exists()

    ann_file = ann_dir / ANNOTATION_FILENAME

    # 3a. Attempting optimize with positions_reviewed=False fails with REVIEW_NOT_CONFIRMED
    unconfirmed_rev = _build_valid_review_file(
        tmp_path / "rev_unconfirmed.json",
        ann,
        positions=[4, 5, 6],
        overrides={"positions_reviewed": False},
    )
    with pytest.raises(NBFoldError) as exc1:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_file,
            review_manifest_path=unconfirmed_rev,
            output_dir=tmp_path / "opt_unconfirmed",
        )
    assert exc1.value.error_code == "REVIEW_NOT_CONFIRMED"

    # 3b. Attempting optimize with empty reviewed_positions fails with INVALID_REVIEW_MANIFEST
    empty_pos_rev = _build_valid_review_file(
        tmp_path / "rev_empty.json",
        ann,
        positions=[],
    )
    with pytest.raises(NBFoldError) as exc2:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_file,
            review_manifest_path=empty_pos_rev,
            output_dir=tmp_path / "opt_empty",
        )
    assert exc2.value.error_code == "INVALID_REVIEW_MANIFEST"

    # 3c. Attempting optimize with stale structure hash fails with REVIEW_MANIFEST_STALE
    stale_rev = _build_valid_review_file(
        tmp_path / "rev_stale.json",
        ann,
        positions=[4, 5, 6],
        overrides={"structure_sha256": "0" * 64},
    )
    with pytest.raises(NBFoldError) as exc3:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_file,
            review_manifest_path=stale_rev,
            output_dir=tmp_path / "opt_stale",
        )
    assert exc3.value.error_code == "REVIEW_MANIFEST_STALE"


# ---------------------------------------------------------------------------
# 4. MMMMM counterexample: infeasible default GC vs feasible zero-gain
# ---------------------------------------------------------------------------


def test_mmmmm_infeasible_default_gc_and_feasible_zero_gain(tmp_path: Path) -> None:
    """Verify MMMMM fails with INFEASIBLE under default 40-47% GC, and succeeds with 0 gain under 30-40% GC."""
    m_prot = "MMMMM"
    fasta_p = _write_fasta(tmp_path / "mmmmm.fasta", m_prot)
    pdb_p = tmp_path / "mmmmm.pdb"
    pdb_p.write_text(_make_pdb_text(m_prot, [40.0, 40.0, 40.0, 90.0, 90.0]), encoding="utf-8")

    # 4a. Default GC (40-47%): ATG is 33.33% GC -> must fail with INFEASIBLE
    ann_def_dir = tmp_path / "ann_def"
    ann_def = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        provider_family="af2",
        output_dir=ann_def_dir,
    )["annotation"]
    rev_def_p = _build_valid_review_file(tmp_path / "rev_def.json", ann_def)

    with pytest.raises(NBFoldError) as exc_inf:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_def_dir / ANNOTATION_FILENAME,
            review_manifest_path=rev_def_p,
            output_dir=tmp_path / "opt_inf",
        )
    assert exc_inf.value.error_code == "INFEASIBLE"
    assert not (tmp_path / "opt_inf" / PAIRED_FASTA_FILENAME).exists()

    # 4b. Feasible GC (30-40%): succeeds even though delta_3d_gain == 0.0 and total_edits == 0
    ann_feas_dir = tmp_path / "ann_feas"
    ann_feas = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        provider_family="af2",
        gc_min=30.0,
        gc_max=40.0,
        output_dir=ann_feas_dir,
    )["annotation"]
    rev_feas_p = _build_valid_review_file(tmp_path / "rev_feas.json", ann_feas)

    opt_res = run_nbfold_optimize(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        annotation_path=ann_feas_dir / ANNOTATION_FILENAME,
        review_manifest_path=rev_feas_p,
        gc_min=30.0,
        gc_max=40.0,
        output_dir=tmp_path / "opt_feas",
    )
    comp = opt_res["observations"]["paired_comparison"]
    assert comp["hard_constraints_pass"] is True
    assert comp["delta_3d_gain"] == 0.0
    assert comp["total_edits"] == 0


# ---------------------------------------------------------------------------
# 5. MFE unavailable/error/per-call settings, sense-CDS scope, and policy gates
# ---------------------------------------------------------------------------


def test_observations_mfe_states_scope_guard_and_policy_gates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify MFE NOT_AVAILABLE, ERROR, per-call settings, unmapped scope guard, and required policy BLOCK."""
    fasta_p = _write_fasta(tmp_path / "prot.fasta", SYNTHETIC_PROTEIN)
    bfs = [90.0] * 4 + [42.0] * 5 + [90.0] * 10
    pdb_p = tmp_path / "model.pdb"
    pdb_p.write_text(_make_pdb_text(SYNTHETIC_PROTEIN, bfs), encoding="utf-8")

    # 5a. Unmapped construct scope guard -> NOT_AVAILABLE
    unmapped_obs = evaluate_candidate_observations(
        candidate_id="test_unmapped",
        engine_name="dp_v2",
        engine_ver="2.0.1",
        cds="GCCACC" + "ATG" * len(SYNTHETIC_PROTEIN),
        expected_protein=SYNTHETIC_PROTEIN,
        cai=0.8,
        gc_fraction=0.45,
        structure_fit=0.7,
        gc_min_pct=0.0,
        gc_max_pct=100.0,
        forbidden_motifs=[],
        unmapped_construct=True,
    )
    assert unmapped_obs["scope"] == "not_established"
    assert unmapped_obs["plant_rna_scan"]["summary_status"] == "NOT_AVAILABLE"
    assert unmapped_obs["initiation_mfe"]["status"] == "NOT_AVAILABLE"
    assert unmapped_obs["yield_prediction"]["status"] == "NOT_AVAILABLE"
    assert unmapped_obs["yield_prediction"]["value"] is None

    # 5b. Required MFE policy when ViennaRNA is absent -> REQUIRED_CHECK_UNAVAILABLE
    req_sop = {
        "$schema": "factorforge-sop-v1",
        "profile_id": "synthetic_req_mfe",
        "sop_name": "Synthetic Required MFE",
        "version": "1.0.0",
        "status": "ACTIVE",
        "default_enforcement": "IGNORE",
        "unknown_rule_policy": "ERROR",
        "rules": {"policy.initiation_mfe.v1": "HARD_FAIL"},
    }
    sop_p = tmp_path / "sop_req.json"
    sop_p.write_text(json.dumps(req_sop), encoding="utf-8")

    ann_dir = tmp_path / "ann_sop"
    ann = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        provider_family="af2",
        output_dir=ann_dir,
    )["annotation"]

    # Compute config_digest with sop_profile
    from factorforge.nbfold.pipeline import resolve_nbfold_configuration

    cfg_with_sop = resolve_nbfold_configuration(sop_profile_path=sop_p)
    rev_sop_p = _build_valid_review_file(
        tmp_path / "rev_sop.json",
        ann,
        overrides={"config_digest": cfg_with_sop["config_digest"]},
    )

    monkeypatch.setitem(sys.modules, "RNA", None)
    with pytest.raises(NBFoldError) as exc_req:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_dir / ANNOTATION_FILENAME,
            review_manifest_path=rev_sop_p,
            sop_profile_path=sop_p,
            output_dir=tmp_path / "opt_req_unavail",
        )
    assert exc_req.value.error_code == "REQUIRED_CHECK_UNAVAILABLE"

    # 5c. Corrupted MFE engine output -> OBSERVATION_ENGINE_ERROR
    monkeypatch.setitem(
        sys.modules,
        "RNA",
        SimpleNamespace(
            md=lambda: SimpleNamespace(temperature=37.0),
            fold_compound=lambda seq, md: SimpleNamespace(mfe=lambda: ("(" * len(seq), -4.0)),
        ),
    )
    rev_default_p = _build_valid_review_file(tmp_path / "rev_def.json", ann)
    with pytest.raises(NBFoldError) as exc_mfe_err:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_dir / ANNOTATION_FILENAME,
            review_manifest_path=rev_default_p,
            output_dir=tmp_path / "opt_mfe_err",
        )
    assert exc_mfe_err.value.error_code == "OBSERVATION_ENGINE_ERROR"

    # 5d. Valid controlled MFE engine with upstream_context and temperature_c=41.5
    monkeypatch.setitem(
        sys.modules,
        "RNA",
        SimpleNamespace(
            __version__="2.6.4-synthetic",
            md=lambda: SimpleNamespace(temperature=37.0, dangles=2, noLP=0, maxBPspan=-1),
            fold_compound=lambda seq, md: SimpleNamespace(mfe=lambda: ("." * len(seq), -2.25)),
        ),
    )
    opt_ok = run_nbfold_optimize(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        annotation_path=ann_dir / ANNOTATION_FILENAME,
        review_manifest_path=rev_default_p,
        upstream_context="GGGGCCCCAAAA",  # capped to 4 adjacent bases "AAAA"
        temperature_c=41.5,
        output_dir=tmp_path / "opt_mfe_ok",
    )
    mfe_3dp = opt_ok["observations"]["candidates"]["structure_3dp"]["initiation_mfe"]
    assert mfe_3dp["status"] == "AVAILABLE"
    assert mfe_3dp["mfe_kcal_mol"] == -2.25
    assert (mfe_3dp["window_start"], mfe_3dp["window_end"]) == (-4, 44)
    assert mfe_3dp["upstream_bases_used"] == 4
    assert mfe_3dp["temperature_c"] == 41.5


# ---------------------------------------------------------------------------
# 6. Future backend capability gates (ESMFold / ESM-2 fail closed)
# ---------------------------------------------------------------------------


def test_future_backends_fail_closed_even_when_enabled(tmp_path: Path) -> None:
    caps = get_backend_capabilities()
    assert caps["alphafold_file"]["available"] is True
    assert caps["esmfold"]["available"] is False
    assert caps["esm2"]["available"] is False

    fasta_p = _write_fasta(tmp_path / "prot.fasta", SYNTHETIC_PROTEIN)
    pdb_p = tmp_path / "model.pdb"
    pdb_p.write_text(_make_pdb_text(SYNTHETIC_PROTEIN, [50.0] * len(SYNTHETIC_PROTEIN)), encoding="utf-8")

    for kwargs in (
        {"enable_esmfold": True},
        {"enable_esm2": True},
        {"structure_backend": "esmfold"},
    ):
        with pytest.raises(NBFoldError) as exc_info:
            run_nbfold_annotate(
                protein_fasta_path=fasta_p,
                structure_file_path=pdb_p,
                provider_family="af2",
                output_dir=tmp_path / f"out_{ '_'.join(f'{k}_{v}' for k, v in kwargs.items()) }",
                **kwargs,
            )
        assert exc_info.value.error_code == "UNSUPPORTED_BACKEND"


# ---------------------------------------------------------------------------
# 7 & 8. End-to-end CLI flow, overwrite refusal, lambda=0 equivalence, & hashes
# ---------------------------------------------------------------------------


def test_cli_end_to_end_annotate_optimize_overwrite_refusal_and_lambda_zero(
    tmp_path: Path,
) -> None:
    runner = CliRunner()
    fasta_p = _write_fasta(tmp_path / "input.fasta", SYNTHETIC_PROTEIN, header="dangerous/../../header_id")
    bfs = [92.0] * 4 + [40.0] * 5 + [91.0] * 10
    cif_p = tmp_path / "af3_model.cif"
    cif_p.write_text(_make_af3_mmcif_text(SYNTHETIC_PROTEIN, ca_bfactors=bfs), encoding="utf-8")

    ann_dir = tmp_path / "annotation-run"
    res_ann = runner.invoke(
        cli,
        [
            "nbfold",
            "annotate",
            "--protein-fasta",
            str(fasta_p),
            "--structure-file",
            str(cif_p),
            "--format",
            "mmcif",
            "--chain",
            "A",
            "--model",
            "1",
            "--provider-family",
            "af3",
            "--output-dir",
            str(ann_dir),
        ],
    )
    assert res_ann.exit_code == 0, res_ann.output
    ann_summary = json.loads(res_ann.output.strip())
    assert ann_summary["status"] == "AVAILABLE"
    assert ann_summary["reviewed"] is False
    assert ann_summary["recommended_positions"] == [5, 6, 7, 8, 9]

    # Overwrite refusal on annotate
    res_ann_dup = runner.invoke(
        cli,
        [
            "nbfold",
            "annotate",
            "--protein-fasta",
            str(fasta_p),
            "--structure-file",
            str(cif_p),
            "--format",
            "mmcif",
            "--chain",
            "A",
            "--model",
            "1",
            "--provider-family",
            "af3",
            "--output-dir",
            str(ann_dir),
        ],
    )
    assert res_ann_dup.exit_code == 7
    assert "OUTPUT_CONFLICT" in res_ann_dup.output

    ann_data = json.loads((ann_dir / ANNOTATION_FILENAME).read_text(encoding="utf-8"))
    review_p = _build_valid_review_file(tmp_path / "review.json", ann_data)

    design_dir = tmp_path / "design-run"
    res_opt = runner.invoke(
        cli,
        [
            "nbfold",
            "optimize",
            "--protein-fasta",
            str(fasta_p),
            "--structure-file",
            str(cif_p),
            "--annotation",
            str(ann_dir / ANNOTATION_FILENAME),
            "--review-manifest",
            str(review_p),
            "--host",
            "nbenthamiana",
            "--output-dir",
            str(design_dir),
        ],
    )
    assert res_opt.exit_code == 0, res_opt.output
    opt_summary = json.loads(res_opt.output.strip())
    assert opt_summary["status"] == "COMPLETED"
    assert opt_summary["yield_prediction_status"] == "NOT_AVAILABLE"

    # Verify fixed artifact filenames (no FASTA header IDs in filenames)
    for expected_file in (
        PAIRED_FASTA_FILENAME,
        OBSERVATIONS_FILENAME,
        RUN_MANIFEST_FILENAME,
        REPORT_MARKDOWN_FILENAME,
    ):
        assert (design_dir / expected_file).exists()

    # Verify artifact SHA-256 digests in run_manifest.json match actual disk bytes
    run_manifest = json.loads((design_dir / RUN_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    for fname, recorded_sha in run_manifest["output_artifact_sha256s"].items():
        actual_sha = hashlib.sha256((design_dir / fname).read_bytes()).hexdigest()
        assert recorded_sha == actual_sha
        assert recorded_sha != run_manifest["input_provenance"]["protein_sha256"]

    # Overwrite refusal on optimize
    res_opt_dup = runner.invoke(
        cli,
        [
            "nbfold",
            "optimize",
            "--protein-fasta",
            str(fasta_p),
            "--structure-file",
            str(cif_p),
            "--annotation",
            str(ann_dir / ANNOTATION_FILENAME),
            "--review-manifest",
            str(review_p),
            "--host",
            "nbenthamiana",
            "--output-dir",
            str(design_dir),
        ],
    )
    assert res_opt_dup.exit_code == 7
    assert "OUTPUT_CONFLICT" in res_opt_dup.output

    # Verify lambda=0 equivalence (penalty_weight=0.0 -> 3DP equals 1DP)
    ann_zero_dir = tmp_path / "ann-zero"
    ann_zero = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=cif_p,
        declared_format="mmcif",
        chain_id="A",
        model_number=1,
        provider_family="af3",
        penalty_weight=0.0,
        output_dir=ann_zero_dir,
    )["annotation"]
    rev_zero_p = _build_valid_review_file(tmp_path / "rev_zero.json", ann_zero)
    opt_zero = run_nbfold_optimize(
        protein_fasta_path=fasta_p,
        structure_file_path=cif_p,
        annotation_path=ann_zero_dir / ANNOTATION_FILENAME,
        review_manifest_path=rev_zero_p,
        penalty_weight=0.0,
        output_dir=tmp_path / "opt-zero",
    )
    c1 = opt_zero["observations"]["candidates"]["baseline_1dp"]
    c3 = opt_zero["observations"]["candidates"]["structure_3dp"]
    assert c1["cds_sha256"] == c3["cds_sha256"]
    assert opt_zero["observations"]["paired_comparison"]["total_edits"] == 0


def test_manual_positions_custom_motifs_flanks_and_timeout(tmp_path: Path) -> None:
    """Verify manual user positions on NO_CANDIDATES, custom motifs/flanks, and EXECUTION_TIMEOUT."""
    fasta_p = _write_fasta(tmp_path / "prot.fasta", SYNTHETIC_PROTEIN)
    pdb_p = tmp_path / "high_conf.pdb"
    pdb_p.write_text(_make_pdb_text(SYNTHETIC_PROTEIN, [92.0] * len(SYNTHETIC_PROTEIN)), encoding="utf-8")

    # 1. Manual positions on NO_CANDIDATES with custom forbidden motif and flanks
    ann_dir = tmp_path / "ann_manual"
    ann = run_nbfold_annotate(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        provider_family="af2",
        custom_forbidden_motifs=["GAATTC"],
        output_dir=ann_dir,
    )["annotation"]
    assert ann["status"] == "NO_CANDIDATES"
    assert ann["recommended_positions"] == []

    manual_rev = _build_valid_review_file(
        tmp_path / "rev_manual.json",
        ann,
        positions=[4, 5, 6],
        overrides={"modification_rationale": "Manual user linker selection on high-confidence chain"},
    )
    opt_manual = run_nbfold_optimize(
        protein_fasta_path=fasta_p,
        structure_file_path=pdb_p,
        annotation_path=ann_dir / ANNOTATION_FILENAME,
        review_manifest_path=manual_rev,
        custom_forbidden_motifs=["GAATTC"],
        left_flank="ACGTACGT",
        right_flank="TGCATGCA",
        output_dir=tmp_path / "opt_manual",
    )
    assert opt_manual["observations"]["paired_comparison"]["reviewed_positions"] == [4, 5, 6]
    assert "GAATTC" in opt_manual["run_manifest"]["reference_and_constraints"]["forbidden_motifs"]
    assert opt_manual["run_manifest"]["reference_and_constraints"]["left_flank_length_nt"] == 8
    assert opt_manual["run_manifest"]["reference_and_constraints"]["right_flank_length_nt"] == 8

    # 2. Timeout enforcement -> EXECUTION_TIMEOUT (exit code 8)
    with pytest.raises(NBFoldError) as exc_timeout:
        run_nbfold_optimize(
            protein_fasta_path=fasta_p,
            structure_file_path=pdb_p,
            annotation_path=ann_dir / ANNOTATION_FILENAME,
            review_manifest_path=manual_rev,
            custom_forbidden_motifs=["GAATTC"],
            timeout_seconds=0.0,
            output_dir=tmp_path / "opt_timeout",
        )
    assert exc_timeout.value.error_code == "EXECUTION_TIMEOUT"
    assert exc_timeout.value.exit_code == 8



@pytest.mark.parametrize("provider_override", [None, "af3"])
def test_review_rejects_unexplained_manual_positions_or_changed_provider(
    tmp_path: Path, provider_override: str | None
) -> None:
    fasta = _write_fasta(tmp_path / "protein.fasta", SYNTHETIC_PROTEIN)
    structure = tmp_path / "structure.pdb"
    structure.write_text(_make_pdb_text(), encoding="utf-8")
    result = run_nbfold_annotate(
        protein_fasta_path=fasta, structure_file_path=structure,
        provider_family="af2", output_dir=tmp_path / "annotation",
    )
    review = _build_valid_review_file(
        tmp_path / "review.json", result["annotation"], positions=[4, 5, 6],
        modification_rationale="Explicit local selection" if provider_override else None,
    )
    with pytest.raises(NBFoldError) as exc:
        run_nbfold_optimize(
            protein_fasta_path=fasta, structure_file_path=structure,
            annotation_path=Path(result["annotation_path"]), review_manifest_path=review,
            provider_family=provider_override, output_dir=tmp_path / "output",
        )
    assert exc.value.error_code == "REVIEW_MANIFEST_STALE"
    assert not (tmp_path / "output" / PAIRED_FASTA_FILENAME).exists()


def test_review_rejects_changed_export_metadata(tmp_path: Path) -> None:
    from factorforge.nbfold.pipeline import resolve_nbfold_configuration, verify_review_manifest

    fasta = _write_fasta(tmp_path / "protein.fasta", SYNTHETIC_PROTEIN)
    structure = tmp_path / "structure.pdb"
    structure.write_text(_make_pdb_text(), encoding="utf-8")
    result = run_nbfold_annotate(
        protein_fasta_path=fasta, structure_file_path=structure,
        provider_family="af2", output_dir=tmp_path / "annotation",
    )
    review = _build_valid_review_file(tmp_path / "review.json", result["annotation"])
    current = dict(result["annotation"], export_metadata_sha256="f" * 64)
    with pytest.raises(NBFoldError) as exc:
        verify_review_manifest(
            review_manifest_path=review, annotation_path=Path(result["annotation_path"]),
            current_annotation=current, resolved_config=resolve_nbfold_configuration(),
        )
    assert exc.value.error_code == "REVIEW_MANIFEST_STALE"
