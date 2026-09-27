# factorforge/tests/test_partial_dp_rescue.py
"""Unit Tests for Local Partial DP Rescue & Hybrid Stitching Engine (Job 286A)."""

from factorforge.discovery.filter import HardConstraintConfig
from factorforge.engines.sllm.compiler import DesignContractCompiler
from factorforge.engines.sllm.decoder import ConstrainedBeamDecoder, ConstrainedSampleDecoder


def test_partial_dp_rescue_stitching_validity():
    contract = DesignContractCompiler.resolve_contract(
        host="nbenthamiana",
        config=HardConstraintConfig(
            forbidden_type_iis_enzymes={"BsaI", "BsmBI"},
            forbidden_other_motifs=["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"],
            homopolymer_max_run=5,
            global_gc_min=30.0,
            global_gc_max=60.0,
            initiation_45nt_gc_min=20.0,
            initiation_45nt_gc_max=50.0,
        ),
        expected_stop="TAA",
    )

    decoder = ConstrainedBeamDecoder(contract=contract)
    protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"

    # Valid initiation prefix
    simulated_prefix = [
        "ATG",
        "TCA",
        "AAA",
        "GGA",
        "GAA",
        "GAA",
        "TTA",
        "TTT",
        "ACA",
        "GGA",
        "GTT",
        "GTT",
        "CCA",
        "ATA",
        "TTA",
        "GTT",
        "GAA",
        "TTA",
        "GAT",
        "GGA",
    ]
    deadend_pos = 20

    res = decoder._execute_partial_dp_rescue(
        protein_sequence=protein,
        stop_codon="TAA",
        deadend_pos=deadend_pos,
        prefix_codons=simulated_prefix,
        lookback_codons=5,
    )

    assert res["status"] == "PARTIAL_DP_RESCUED"
    assert res["rescue_triggered"] is True
    assert res["partial_rescue"] is True
    assert res["preserved_prefix_codons"] == 15  # 20 - 5
    assert res["sequence"].startswith("".join(simulated_prefix[:15]))
    assert res["sequence"].endswith("TAA")

    # Verify 100% feasibility under independent final_validator
    val = contract.final_validator.evaluate(
        sequence_dna=res["sequence"],
        expected_protein_aa=protein,
        expected_stop="TAA",
    )
    assert val.is_feasible is True, f"Violations: {val.violations}"


def test_partial_dp_rescue_invalid_prefix_triggers_global_fallback():
    contract = DesignContractCompiler.resolve_contract(
        host="nbenthamiana",
        config=HardConstraintConfig(
            forbidden_type_iis_enzymes={"BsaI", "BsmBI"},
            homopolymer_max_run=5,
            initiation_45nt_gc_min=15.0,
            initiation_45nt_gc_max=45.0,
        ),
        expected_stop="TAA",
    )

    decoder = ConstrainedBeamDecoder(contract=contract)
    protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"

    # High-GC prefix that violates 5' initiation contract
    bad_gc_prefix = [
        "ATG",
        "AGC",
        "AAG",
        "GGA",
        "GAG",
        "GAG",
        "CTG",
        "TTC",
        "ACT",
        "GGT",
        "GTT",
        "GTC",
        "CCT",
        "ATC",
        "CTG",
        "GTT",
        "GAG",
        "CTG",
        "GAT",
        "GGT",
    ]

    res = decoder._execute_partial_dp_rescue(
        protein_sequence=protein,
        stop_codon="TAA",
        deadend_pos=20,
        prefix_codons=bad_gc_prefix,
        lookback_codons=5,
    )

    # When partial stitch violates the contract, it MUST safely fall back to global DP rescue
    assert res["status"] == "DP_RESCUED"
    assert res["rescue_triggered"] is True
    assert any("validation_failed" in item for item in res["fallback_chain"])

    val = contract.final_validator.evaluate(
        sequence_dna=res["sequence"],
        expected_protein_aa=protein,
        expected_stop="TAA",
    )
    assert val.is_feasible is True


def test_partial_dp_rescue_early_deadend_fallback():
    contract = DesignContractCompiler.resolve_contract(
        host="nbenthamiana",
        config=HardConstraintConfig(
            forbidden_type_iis_enzymes={"BsaI", "BsmBI"},
            homopolymer_max_run=5,
        ),
        expected_stop="TAA",
    )

    decoder = ConstrainedBeamDecoder(contract=contract)
    protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"

    # Dead-end too early (pos=2)
    res = decoder._execute_partial_dp_rescue(
        protein_sequence=protein,
        stop_codon="TAA",
        deadend_pos=2,
        prefix_codons=["ATG", "AGC"],
        lookback_codons=5,
    )

    # Must safely fallback to global DP rescue
    assert res["status"] == "DP_RESCUED"
    assert res["rescue_triggered"] is True

    val = contract.final_validator.evaluate(
        sequence_dna=res["sequence"],
        expected_protein_aa=protein,
        expected_stop="TAA",
    )
    assert val.is_feasible is True


def test_sample_decoder_partial_rescue_integration():
    contract = DesignContractCompiler.resolve_contract(
        host="nbenthamiana",
        config=HardConstraintConfig(
            forbidden_type_iis_enzymes={"BsaI", "BsmBI"},
            forbidden_other_motifs=["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"],
            homopolymer_max_run=5,
            global_gc_min=35.0,
            global_gc_max=55.0,
        ),
        expected_stop="TAA",
    )

    sample_decoder = ConstrainedSampleDecoder(
        contract=contract,
        temperature=0.8,
        seed=123,
    )
    target_aa = "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEAFQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGCFIYKVKFIGVNFPSDGPVMQKKTMGWEASTERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"

    cand = sample_decoder.decode(target_aa, stop_codon="TAA")
    assert cand["status"] in ("SUCCESS", "PARTIAL_DP_RESCUED", "DP_RESCUED")

    val = contract.final_validator.evaluate(
        sequence_dna=cand["sequence"],
        expected_protein_aa=target_aa,
        expected_stop="TAA",
    )
    assert val.is_feasible is True
