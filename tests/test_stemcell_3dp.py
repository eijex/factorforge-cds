"""Tests for FactorForge 3DP-SC: Human Stem Cell Optimization Engine."""

import pytest
from factorforge.engines.stemcell_3dp import (
    StemCell3DPOptimizer,
    load_stemcell_codon_weights,
    evaluate_37c_thermodynamics,
    calculate_max_duplex_stem,
    STEM_CELL_DEFAULT_FORBIDDEN_MOTIFS,
)


def test_stemcell_codon_weights():
    weights = load_stemcell_codon_weights()
    assert len(weights) == 64  # Complete universal 64 codons
    assert all(0 < w <= 1.0 for w in weights.values())
    # Common high-efficiency codons should have high weight
    assert weights["ATG"] == 1.0  # Met
    assert weights["TGG"] == 1.0  # Trp


def test_max_duplex_stem_calculation():
    # 5 bp stem: (((((...)))))
    struct1 = "(((((...)))))"
    assert calculate_max_duplex_stem(struct1) == 5

    # 22 bp stem (PKR risk): 22 opening brackets
    struct2 = "(" * 22 + "...." + ")" * 22
    assert calculate_max_duplex_stem(struct2) == 22


def test_stemcell_egfp_optimization():
    pytest.importorskip("RNA")
    egfp_prot = (
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
        "RYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNY"
        "NSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHM"
        "VLLEFVTAAGITLGMDELYK"
    )

    opt = StemCell3DPOptimizer(target_gc_min=0.44, target_gc_max=0.50)
    res = opt.optimize(egfp_prot)

    assert res["status"] == "REVIEW_REQUIRED"
    assert res["translation_identity"] == 1.0
    assert res["dna_length_nt"] == len(egfp_prot) * 3
    assert 44.0 <= res["gc_percent"] <= 50.0
    assert res["cai"] > 0.90
    assert res["motif_clean"] is True
    assert res["thermodynamics_37c"]["pkr_assessment"]["status"] == "NOT_AVAILABLE"
    assert res["thermodynamics_37c"]["max_duplex_stem_bp"] < 20

    # Verify no miRNA seeds or Type IIS motifs exist in sequence
    seq = res["optimized_sequence"]
    for motif in STEM_CELL_DEFAULT_FORBIDDEN_MOTIFS:
        assert motif not in seq, f"Forbidden motif {motif} leaked into sequence!"


def test_stemcell_cd9_optimization():
    pytest.importorskip("RNA")
    cd9_prot = (
        "MPVKGGTKCIKYLLFGFNFIFWLAGIAVLAIGLWLRFDSQTKSIFEQETNNNNSSFYTGVYILIGAGALMMLV"
        "GFLGCCGAVQESQCMLGLFFGFLLVIFAIEIAAAIWGYSHKDEVIKEVQEFYKDTYNKLKTKDEPQRETLKAI"
        "HYALNCCGLAGGVEQFISDICPKKDVLETFTVKSCPDAIKEVFDNKFHIIGAVGIGIAVVMIFGMIFSMILCC"
        "AIRRNREMV"
    )

    opt = StemCell3DPOptimizer(target_gc_min=0.45, target_gc_max=0.50)
    res = opt.optimize(cd9_prot)

    assert res["status"] == "REVIEW_REQUIRED"
    assert res["translation_identity"] == 1.0
    assert 45.0 <= res["gc_percent"] <= 50.0
    assert res["motif_clean"] is True
    assert res["thermodynamics_37c"]["pkr_assessment"]["status"] == "NOT_AVAILABLE"


@pytest.mark.parametrize(
    "structure,expected",
    [
        ("(" * 20 + "..." + ")." * 20, 1),
        ("(((...)))...((...))", 3),
        ("...", 0),
        ("((..((...))..))", 2),
    ],
)
def test_stacked_pairs_with_bulges_and_branches(structure, expected):
    assert calculate_max_duplex_stem(structure) == expected


@pytest.mark.parametrize("structure", ["(()", ")(", "[]", "(.x.)"])
def test_bad_structures_rejected(structure):
    with pytest.raises(ValueError):
        calculate_max_duplex_stem(structure)


def test_temperature_is_per_call():
    RNA = pytest.importorskip("RNA")
    original = RNA.cvar.temperature
    try:
        RNA.cvar.temperature = 25.0
        result = evaluate_37c_thermodynamics("ATGGCTGCT")
        assert result["status"] == "AVAILABLE"
        assert RNA.cvar.temperature == 25.0
        assert result["temperature_c"] == 37.0
    finally:
        RNA.cvar.temperature = original


def test_missing_engine_never_passes(monkeypatch):
    import sys

    monkeypatch.setitem(sys.modules, "RNA", None)
    result = StemCell3DPOptimizer(target_gc_min=0, target_gc_max=1).optimize("MAAAA")
    assert result["computational_status"] == "INCOMPLETE"
    assert result["checks"]["initiation_mfe"] is None
    assert result["thermodynamics_37c"]["status"] == "NOT_AVAILABLE"


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), True])
def test_invalid_engine_output_is_error(monkeypatch, bad):
    import sys
    from types import SimpleNamespace

    fake = SimpleNamespace(
        md=lambda: SimpleNamespace(temperature=25),
        fold_compound=lambda seq, md: SimpleNamespace(mfe=lambda: ("." * len(seq), bad)),
    )
    monkeypatch.setitem(sys.modules, "RNA", fake)
    assert evaluate_37c_thermodynamics("ATGGCT")["status"] == "ERROR"


def test_thermo_failure_blocks_computational_status(monkeypatch):
    import factorforge.engines.stemcell_3dp as module

    monkeypatch.setattr(
        module,
        "evaluate_37c_thermodynamics",
        lambda dna: dict(
            status="AVAILABLE",
            initiation_mfe_37c=-100.0,
            min_local_mfe_40nt=-100.0,
            max_duplex_stem_bp=1,
        ),
    )
    result = StemCell3DPOptimizer(target_gc_min=0, target_gc_max=1).optimize("MAAAA")
    assert result["computational_status"] == "BLOCKED"
    assert result["checks"]["initiation_mfe"] is False
    assert result["checks"]["local_mfe"] is False
    assert result["relaxation_performed"] is False


def test_review_required_for_structure_positions():
    with pytest.raises(ValueError, match="reviewed"):
        StemCell3DPOptimizer(pause_positions=[2])


def test_zero_penalty_matches_same_reference_baseline():
    base = StemCell3DPOptimizer(target_gc_min=0, target_gc_max=1)
    other = StemCell3DPOptimizer(
        target_gc_min=0,
        target_gc_max=1,
        pause_positions=[2, 3],
        positions_reviewed=True,
        structure_reference="synthetic_reference",
        pause_penalty_weight=0,
    )
    assert (
        base.optimize("MASTWQHLEK")["optimized_sequence"]
        == other.optimize("MASTWQHLEK")["optimized_sequence"]
    )


def test_policy_is_registered_packaged_and_consumed():
    import json
    from pathlib import Path
    from factorforge.engines.stemcell_3dp import POLICY_PATH
    from factorforge.rules.registry import RuleRegistry
    from factorforge.sops.config import SopProfile

    config = Path(__file__).parents[1] / "configs/sops/human_stemcell_conservative_reference.json"
    profile = SopProfile(config, known_rule_ids=[r.rule_id for r in RuleRegistry().list_rules()])
    assert json.loads(config.read_text()) == json.loads(POLICY_PATH.read_text())
    engine = StemCell3DPOptimizer(target_gc_min=0, target_gc_max=1)
    assert engine.profile.digest == profile.digest
    result = engine.evaluate_candidate("ATGGCTAATAAAGCT", "MANKA")
    assert result["computational_status"] == "BLOCKED"
    assert result["registered_policy"]["hard_fail_count"] >= 1


@pytest.mark.parametrize(
    "minimum,maximum", [(float("nan"), 0.5), (0.6, 0.4), (-1, 0.5), (False, 0.5)]
)
def test_invalid_gc_bounds_rejected(minimum, maximum):
    with pytest.raises(ValueError):
        StemCell3DPOptimizer(target_gc_min=minimum, target_gc_max=maximum)


def test_reference_is_explicitly_uncalibrated():
    import json
    from pathlib import Path

    path = Path(__file__).parents[1] / "src/factorforge/data/human_stemcell_codons.json"
    data = json.loads(path.read_text())
    assert data["calibration_status"] == "UNVERIFIED"
    assert data["laboratory_approval"] == "NOT_OBTAINED"


@pytest.mark.parametrize("structure", ["(", ")", "x", "...("])
def test_fold_rejects_malformed_structure(monkeypatch, structure):
    import sys
    from types import SimpleNamespace

    def fold(text, md):
        return SimpleNamespace(mfe=lambda: (structure.ljust(len(text), ".")[: len(text)], -1.0))

    monkeypatch.setitem(
        sys.modules,
        "RNA",
        SimpleNamespace(md=lambda: SimpleNamespace(temperature=25), fold_compound=fold),
    )
    assert evaluate_37c_thermodynamics("ATGGCT")["status"] == "ERROR"


def test_independent_gc_check_blocks_candidate():
    result = StemCell3DPOptimizer().evaluate_candidate("ATGGCTGCT", "MAA")
    assert result["checks"]["gc_in_band"] is False
    assert result["computational_status"] == "BLOCKED"


def test_same_reference_1dp_comparator():
    from factorforge.engines.dp_v2_1_1 import DPV211Optimizer

    engine = StemCell3DPOptimizer(target_gc_min=0, target_gc_max=1)
    expected = DPV211Optimizer(
        forbidden_motifs=engine.forbidden_motifs, include_reverse_complement=False
    ).optimize(
        protein_sequence="MASTWQHLEK",
        codon_weights=engine.codon_weights,
        target_gc_min=0,
        target_gc_max=1,
    )
    assert engine.optimize("MASTWQHLEK")["optimized_sequence"] == expected["sequence"]


def test_profile_thresholds_are_validated(tmp_path):
    import json
    from factorforge.engines.stemcell_3dp import POLICY_PATH

    data = json.loads(POLICY_PATH.read_text())
    data["experimental_observation_policy"]["initiation_min_kcal_mol"] = "not_a_number"
    path = tmp_path / "bad_policy.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="thresholds"):
        StemCell3DPOptimizer(policy_path=path)
