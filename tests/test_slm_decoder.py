# factorforge/tests/test_slm_decoder.py
import pytest
import numpy as np
import math

from factorforge.analysis.metrics import translate_dna
from factorforge.discovery.filter import HardConstraintConfig, HardConstraintFilter
from factorforge.engines.sllm.adapters import DeterministicMockSLMAdapter, SLMModelAdapter
from factorforge.engines.sllm.compiler import DesignContractCompiler
from factorforge.engines.sllm.decoder import ConstrainedBeamDecoder, ConstrainedSampleDecoder
from factorforge.engines.sllm.interfaces import CodonLogits


class AdversarialBlockedAdapter(SLMModelAdapter):
    """Adapter that intentionally returns -inf for all valid codons at step 2 to force dead-end."""

    def __init__(self):
        pass

    @property
    def model_name(self) -> str:
        return "adversarial_blocked_adapter"

    def predict_logits(self, prefix_codons, target_aa, position) -> CodonLogits:
        scores = np.zeros(64, dtype=np.float32)
        if position == 2:
            # Force everything to -inf at pos 2
            scores[:] = -math.inf
        return CodonLogits(scores)


def test_design_contract_compiler():
    config = HardConstraintConfig(
        homopolymer_max_run=4,
        forbidden_other_motifs=["GCGGCCGC", "GGTCTC"],
    )
    pipeline = DesignContractCompiler.compile_from_config(config)
    assert len(pipeline.processors) == 3


def test_constrained_beam_decoder_exact_translation():
    target_protein = "MDLSTLSQPVVSVEP"
    config = HardConstraintConfig()
    filter_guard = HardConstraintFilter(config)

    decoder = ConstrainedBeamDecoder(
        model_adapter=DeterministicMockSLMAdapter(host="nbenthamiana", seed=42),
        beam_width=8,
        hard_config=config,
    )

    res = decoder.decode(target_protein, stop_codon="TAA")
    assert res["candidate_id"] == "cand_sllm_beam_best"

    seq = res["sequence"]
    # Check 100% exact translation
    translated = translate_dna(seq)
    assert translated == target_protein + "*"

    # Run through HardConstraintFilter
    filter_res = filter_guard.evaluate(seq, target_protein)
    assert filter_res.is_feasible
    assert len(filter_res.violations) == 0


def test_constrained_beam_decoder_bsai_veto():
    # Target containing Glycine and Leucine (which easily form GGTCTC BsaI if unconstrained)
    target_protein = "MGLGLGLGL"
    config = HardConstraintConfig(forbidden_other_motifs=["GGTCTC"])
    filter_guard = HardConstraintFilter(config)

    decoder = ConstrainedBeamDecoder(
        beam_width=8,
        hard_config=config,
    )
    res = decoder.decode(target_protein, stop_codon="TAA")
    seq = res["sequence"]

    assert "GGTCTC" not in seq
    assert "GAGACC" not in seq  # Reverse complement

    filter_res = filter_guard.evaluate(seq, target_protein)
    assert filter_res.is_feasible


def test_constrained_sample_decoder_diversity():
    target_protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    config = HardConstraintConfig()

    sample_decoder1 = ConstrainedSampleDecoder(temperature=0.7, hard_config=config, seed=10)
    sample_decoder2 = ConstrainedSampleDecoder(temperature=0.7, hard_config=config, seed=99)

    res1 = sample_decoder1.decode(target_protein[:25], stop_codon="TAA")
    res2 = sample_decoder2.decode(target_protein[:25], stop_codon="TAA")

    assert res1["status"] == "SUCCESS"
    assert res2["status"] == "SUCCESS"
    # Both translate accurately
    assert translate_dna(res1["sequence"]) == target_protein[:25] + "*"
    assert translate_dna(res2["sequence"]) == target_protein[:25] + "*"


def test_exact_dp_rescue_on_forced_dead_end():
    target_protein = "MDLSTLSQPVVSVEP"
    config = HardConstraintConfig()

    adversarial_adapter = AdversarialBlockedAdapter()
    decoder = ConstrainedBeamDecoder(
        model_adapter=adversarial_adapter,
        beam_width=2,
        hard_config=config,
    )

    res = decoder.decode(target_protein, stop_codon="TAA")
    # Must have triggered exact DP rescue
    assert res["status"] == "DP_RESCUED"
    assert res["rescue_triggered"]
    assert res["candidate_id"] == "cand_sllm_dp_rescue"
    assert res["requested_generator"] == "sllm_hybrid_beam"
    assert res["actual_generator"] == "dp_v2_1_1"
    assert res["score"] is None
    assert res["score_type"] == "dp_objective"
    assert "dp_v2_1_1_rescue" in res["fallback_chain"]
    assert "Exact DP rescue executed" in res["rationale"]

    # Check that rescued sequence is 100% valid
    filter_guard = HardConstraintFilter(config)
    filter_res = filter_guard.evaluate(res["sequence"], target_protein, expected_stop="TAA")
    assert filter_res.is_feasible


def test_expected_stop_mismatch_fails_hard_filter():
    filter_guard = HardConstraintFilter()
    target_protein = "MKLL"
    valid_seq_taa = "ATGAAGCTGCTCTAA"
    valid_seq_tag = "ATGAAGCTGCTCTAG"
    no_stop_seq = "ATGAAGCTGCTC"

    # Test matching expected stop
    assert filter_guard.evaluate(valid_seq_taa, target_protein, expected_stop="TAA").is_feasible

    # Test mismatching stop codon
    res_mismatch = filter_guard.evaluate(valid_seq_tag, target_protein, expected_stop="TAA")
    assert not res_mismatch.is_feasible
    assert any("expected stop codon 'TAA'" in v for v in res_mismatch.violations)

    # Test missing stop codon
    res_missing = filter_guard.evaluate(no_stop_seq, target_protein, expected_stop="TAA")
    assert not res_missing.is_feasible


def test_internal_stop_codon_fails_hard_filter():
    filter_guard = HardConstraintFilter()
    target_protein = "MKLL"
    internal_stop_seq = "ATGAAGTAACTCTAA"  # TAA at codon 3 is internal stop

    res = filter_guard.evaluate(internal_stop_seq, target_protein, expected_stop="TAA")
    assert not res.is_feasible
    assert any("Internal stop codons detected" in v for v in res.violations)


def test_flank_junction_forbidden_motif_detected():
    config = HardConstraintConfig(forbidden_other_motifs=["GGTCTC"])
    filter_guard = HardConstraintFilter(config)

    # Upstream ends with GG, CDS starts with TCTC -> Forms GGTCTC across junction
    upstream = "AAAAAGG"
    cds = "TCTCATGAAGCTGCTCTAA"  # (Hypothetical)

    res = filter_guard.evaluate(
        sequence_dna=cds,
        expected_protein_aa="SHEALL",
        expected_stop="TAA",
        upstream_flank=upstream,
    )
    assert not res.is_feasible
    assert any("Forbidden restriction motifs found" in v for v in res.violations)


def test_dp_rescue_strict_infeasible_raises_unsat():
    # Set impossible GC limits (e.g. 5% - 10%) so DP cannot find feasible sequence
    impossible_config = HardConstraintConfig(
        global_gc_min=5.0,
        global_gc_max=10.0,
    )
    adversarial_adapter = AdversarialBlockedAdapter()
    decoder = ConstrainedBeamDecoder(
        model_adapter=adversarial_adapter,
        beam_width=2,
        hard_config=impossible_config,
    )

    with pytest.raises(RuntimeError) as exc_info:
        decoder.decode("MDLSTLSQPVVSVEP", stop_codon="TAA")

    assert "UNSAT" in str(exc_info.value)


def test_sample_decoder_fallback_chain_provenance():
    adversarial_adapter = AdversarialBlockedAdapter()
    sample_decoder = ConstrainedSampleDecoder(
        model_adapter=adversarial_adapter,
    )
    res = sample_decoder.decode("MDLSTLSQPVVSVEP", stop_codon="TAA")

    assert res["status"] == "DP_RESCUED"
    assert res["requested_generator"] == "sllm_hybrid_sample"
    assert res["actual_generator"] == "dp_v2_1_1"
    assert "sllm_sample_dead_end" in res["fallback_chain"]
    assert "dp_v2_1_1_rescue" in res["fallback_chain"]


def test_beam_search_normalized_log_prob_score():
    decoder = ConstrainedBeamDecoder()
    res = decoder.decode("MDLSTLSQPVVSVEP", stop_codon="TAA")

    assert res["status"] == "SUCCESS"
    assert res["score_type"] == "constrained_sequence_log_likelihood"
    # Log probabilities are non-positive
    assert res["score"] <= 0.0
