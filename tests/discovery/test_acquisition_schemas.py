"""Unit tests for Paired DBTL data models and serialization (Job 283C)."""

from factorforge.discovery.acquisition import (
    ConstructDesignRecord,
    DerivedOutcomeRecord,
    ExperimentRunRecord,
    MeasurementRecord,
    PairedDBTLDataset,
    SampleRecord,
)
from factorforge.discovery.schemas import TraitVector


def test_construct_design_record_serialization():
    tv = TraitVector(
        cai_golden_set=0.88,
        global_gc_percent=45.2,
        initiation_mfe_kcal_mol=-5.5,
        local_50bp_gc_min=38.0,
        local_50bp_gc_max=52.0,
        homopolymer_max_run=3,
    )
    rec = ConstructDesignRecord(
        construct_id="sfGFP_H0",
        target_name="sfGFP",
        mature_protein_aa_length=238,
        construct_aa_length=238,
        hypothesis_id="H0",
        hypothesis_name="Primary Deterministic Optimum",
        generation_contract="dp_v2_1_1_deterministic_optimum",
        sequence_digest="abc123sha256",
        trait_vector=tv,
        sequence_dna="ATGGCC...",
        is_control=False,
    )

    d_full = rec.to_dict(include_sequence=True)
    assert d_full["sequence_dna"] == "ATGGCC..."
    assert d_full["trait_vector"]["cai_golden_set"] == 0.88

    d_masked = rec.to_dict(include_sequence=False)
    assert "sequence_dna" not in d_masked
    assert d_masked["sequence_digest"] == "abc123sha256"


def test_paired_dbtl_dataset_roundtrip(tmp_path):
    tv = TraitVector(
        cai_golden_set=0.85,
        global_gc_percent=44.0,
        initiation_mfe_kcal_mol=-4.2,
        local_50bp_gc_min=39.0,
        local_50bp_gc_max=50.0,
        homopolymer_max_run=3,
    )
    c_rec = ConstructDesignRecord(
        construct_id="sfGFP_H0",
        target_name="sfGFP",
        mature_protein_aa_length=238,
        construct_aa_length=238,
        hypothesis_id="H0",
        hypothesis_name="Primary Deterministic Optimum",
        generation_contract="dp_v2_1_1_deterministic_optimum",
        sequence_digest="sha256_mock_digest",
        trait_vector=tv,
        sequence_dna="ATGTCTAAAGGT...",
    )
    e_rec = ExperimentRunRecord(
        experiment_id="EXP-001",
        batch_id="BATCH-01",
    )
    s_rec = SampleRecord(
        sample_id="SMP-001",
        blinded_plate_position="A01",
        experiment_id="EXP-001",
        construct_id="sfGFP_H0",
        biological_replicate=1,
    )
    m_rec = MeasurementRecord(
        measurement_id="MSR-001",
        sample_id="SMP-001",
        assay_type="GFP_FLUORESCENCE",
        raw_yield_ug_g_fw=142.5,
    )
    o_rec = DerivedOutcomeRecord(
        outcome_id="OUT-001",
        sample_id="SMP-001",
        construct_id="sfGFP_H0",
        experiment_id="EXP-001",
        normalized_yield_ug_g_fw=142.5,
        relative_to_h0=1.0,
        relative_to_pos_ctrl=1.15,
        lod_loq_status="ABOVE_LOQ",
    )

    dataset = PairedDBTLDataset(
        dataset_id="DBTL-EXP-001",
        constructs={"sfGFP_H0": c_rec},
        experiments={"EXP-001": e_rec},
        samples={"SMP-001": s_rec},
        measurements=[m_rec],
        derived_outcomes=[o_rec],
    )
    dataset.archive_sha256 = dataset.compute_sha256()

    data_dict = dataset.to_dict(include_sequence=True)
    reconstituted = PairedDBTLDataset.from_dict(data_dict)

    assert reconstituted.dataset_id == "DBTL-EXP-001"
    assert "sfGFP_H0" in reconstituted.constructs
    assert reconstituted.constructs["sfGFP_H0"].trait_vector.cai_golden_set == 0.85
    assert len(reconstituted.measurements) == 1
    assert reconstituted.measurements[0].raw_yield_ug_g_fw == 142.5
    assert len(reconstituted.derived_outcomes) == 1
    assert reconstituted.derived_outcomes[0].relative_to_pos_ctrl == 1.15
