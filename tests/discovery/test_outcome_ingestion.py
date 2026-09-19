"""Unit tests for wet-lab outcome ingestion and derived metrics calculation (Job 283C)."""

from factorforge.discovery.acquisition_logger import AcquisitionLogger
from factorforge.discovery.panel_builder import ProspectivePanelBuilder


def test_outcome_ingestion_and_fold_change_calculation():
    builder = ProspectivePanelBuilder(host="nbenthamiana")
    dataset, aux_data = builder.build_panel(
        replicates_per_construct=3,
        seed=42,
        experiment_id="EXP-TEST-INGEST-01",
    )

    # Prepare mock measurements for all 33 samples
    # Let sfGFP_H0 average 120 µg/g, sfGFP_H1 average 150 µg/g, POS_CTRL_sfGFP average 100 µg/g, NEG_CTRL average 0.1 µg/g
    mock_measurements = []
    unblinded = {m["sample_id"]: m for m in aux_data["unblinded_mapping"]}

    for idx, (smp_id, smp_rec) in enumerate(dataset.samples.items()):
        info = unblinded[smp_id]
        cid = info["construct_id"]
        rep = info["biological_replicate"]

        if cid == "sfGFP_H0":
            raw_yield = 120.0 + (rep * 2.0)
        elif cid == "sfGFP_H1":
            raw_yield = 150.0 + (rep * 3.0)
        elif cid == "sfGFP_H2":
            raw_yield = 110.0 - (rep * 1.5)
        elif cid == "POS_CTRL_sfGFP":
            raw_yield = 100.0 + (rep * 1.0)
        elif cid == "NEG_CTRL_EMPTY":
            raw_yield = 0.1
        else:
            raw_yield = 45.0 + rep

        mock_measurements.append(
            {
                "measurement_id": f"MSR-{idx + 1:03d}",
                "sample_id": smp_id,
                "assay_type": "GFP_FLUORESCENCE" if "sfGFP" in cid else "QUANTITATIVE_ELISA",
                "raw_yield_ug_g_fw": raw_yield,
                "qc_status": "PASS",
            }
        )

    # Ingest measurements
    updated_dataset = AcquisitionLogger.ingest_wet_lab_measurements(
        dataset=dataset,
        raw_measurements=mock_measurements,
    )

    assert len(updated_dataset.measurements) == 33
    assert len(updated_dataset.derived_outcomes) == 33

    # Check sfGFP_H1 outcomes: should have relative_to_h0 > 1.0 and relative_to_pos_ctrl > 1.0
    h1_outcomes = [o for o in updated_dataset.derived_outcomes if o.construct_id == "sfGFP_H1"]
    assert len(h1_outcomes) == 3
    for o in h1_outcomes:
        assert o.relative_to_h0 > 1.0
        assert o.relative_to_pos_ctrl > 1.0
        assert o.lod_loq_status == "ABOVE_LOQ"

    # Check NEG_CTRL: should be below LOD (< 0.5)
    neg_outcomes = [
        o for o in updated_dataset.derived_outcomes if o.construct_id == "NEG_CTRL_EMPTY"
    ]
    assert len(neg_outcomes) == 3
    for o in neg_outcomes:
        assert o.lod_loq_status == "BELOW_LOD"
