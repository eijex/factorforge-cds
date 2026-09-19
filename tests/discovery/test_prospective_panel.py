"""Unit tests for ProspectivePanelBuilder, OrthogonalityGate, and Package Export (Job 283C)."""

import json

from factorforge.discovery.acquisition_logger import AcquisitionLogger
from factorforge.discovery.panel_builder import (
    OrthogonalityGate,
    ProspectivePanelBuilder,
)


def test_orthogonality_gate_calculation():
    from factorforge.discovery.schemas import TraitVector

    tv1 = TraitVector(
        cai_golden_set=0.92,
        global_gc_percent=48.0,
        initiation_mfe_kcal_mol=-8.0,
        local_50bp_gc_min=40.0,
        local_50bp_gc_max=55.0,
        homopolymer_max_run=3,
    )
    tv2 = TraitVector(
        cai_golden_set=0.75,
        global_gc_percent=38.0,
        initiation_mfe_kcal_mol=-2.0,
        local_50bp_gc_min=32.0,
        local_50bp_gc_max=45.0,
        homopolymer_max_run=3,
    )

    gate = OrthogonalityGate(min_epsilon=0.05)
    d = gate.compute_distance(tv1, tv2)
    assert d > 0.05

    traits = {"H0": tv1, "H1": tv2}
    passed, dists = gate.verify_orthogonality(traits)
    assert passed is True
    assert "H0_vs_H1" in dists


def test_prospective_panel_builder_execution(tmp_path):
    builder = ProspectivePanelBuilder(host="nbenthamiana")

    # Test on standard 3-target panel
    dataset, aux_data = builder.build_panel(
        replicates_per_construct=3,
        seed=42,
        experiment_id="EXP-TEST-001",
    )

    # 3 targets x 3 hypotheses = 9 constructs + 2 controls = 11 constructs
    assert len(dataset.constructs) == 11
    assert "sfGFP_H0" in dataset.constructs
    assert "sfGFP_H1" in dataset.constructs
    assert "sfGFP_H2" in dataset.constructs
    assert "VP28_H0" in dataset.constructs
    assert "CD47_Ectodomain_H0" in dataset.constructs
    assert "POS_CTRL_sfGFP" in dataset.constructs
    assert "NEG_CTRL_EMPTY" in dataset.constructs

    # Verify 11 constructs x 3 replicates = 33 sample wells
    assert len(dataset.samples) == 33
    assert len(aux_data["blinded_layout"]) == 33
    assert len(aux_data["unblinded_mapping"]) == 33

    # Verify Orthogonality Gate passed for all 3 targets
    ortho_rep = aux_data["orthogonality_report"]
    assert "sfGFP" in ortho_rep
    assert ortho_rep["sfGFP"]["orthogonality_passed"] is True
    assert "VP28" in ortho_rep
    assert ortho_rep["VP28"]["orthogonality_passed"] is True
    assert "CD47_Ectodomain" in ortho_rep
    assert ortho_rep["CD47_Ectodomain"]["orthogonality_passed"] is True

    # Test export package
    logger = AcquisitionLogger(base_output_dir=tmp_path)
    pkg_dir = logger.export_prospective_panel(dataset, aux_data, subfolder="test_panel")

    assert (pkg_dir / "panel_sequences.fasta").exists()
    assert (pkg_dir / "synthesis_manifest.csv").exists()
    assert (pkg_dir / "blinded_plate_layout.json").exists()
    assert (pkg_dir / "unblinded_mapping.json").exists()
    assert (pkg_dir / "scientific_memory_panel.json").exists()
    assert (pkg_dir / "prospective_dbtl_dataset.json").exists()
    assert (pkg_dir / "secure_archive_manifest.json").exists()

    # Verify sequence-free memory does not leak sequence_dna
    with (pkg_dir / "scientific_memory_panel.json").open("r", encoding="utf-8") as f:
        mem = json.load(f)
        for cid, cdata in mem["constructs"].items():
            assert "sequence_dna" not in cdata
            assert "sequence_digest" in cdata
