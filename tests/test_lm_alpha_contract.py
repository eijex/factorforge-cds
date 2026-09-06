"""PyTest contract tests for FactorForge-LM (v3.5.0 ML Track).

Verifies control token encoding, Type IIS logit masking, autoregressive
candidate generation, and design package contract adherence.
"""

from factorforge.engines.lm.tokenizer.control_tokenizer import FactorForgeControlTokenizer
from factorforge.engines.lm.adapter import FactorForgeLogitMasker
from factorforge.engines.lm.inference import ConstrainedBeamSearchEngine
from factorforge.analysis.metrics import translate_dna
from factorforge.evaluation.models import EvaluationResult


def test_control_tokenizer_vocab_and_encoding():
    tokenizer = FactorForgeControlTokenizer()
    assert tokenizer.vocab_size == 98

    encoded = tokenizer.encode_encoder_input("MAKW", host="nbenthamiana", gc_band="40-47")
    decoded = [tokenizer.id_to_token[i] for i in encoded]

    assert decoded[0] == "<s>"
    assert decoded[1] == "<host:nbenthamiana>"
    assert decoded[2] == "<gc:40-47>"
    assert decoded[3] == "<type2is:clean>"
    assert decoded[-1] == "</s>"


def test_logit_masker_bsai_exclusion():
    tokenizer = FactorForgeControlTokenizer()
    masker = FactorForgeLogitMasker(tokenizer)

    # Prefix ending in GGTCT
    prefix = "ATGGGTCT"
    sample_logits = [0.0] * tokenizer.vocab_size

    # Codon CGA starts with C and completes BsaI (GGTCT + C -> GGTCTC)
    # Codon AAG starts with A and does NOT complete BsaI (GGTCT + A -> GGTCTA)
    cga_id = tokenizer.token_to_id["CGA"]
    aag_id = tokenizer.token_to_id["AAG"]

    masked = masker.apply_logit_masks(prefix, sample_logits)
    assert masked[cga_id] == -float("inf"), "BsaI restriction site codon was not masked to -inf!"
    assert masked[aag_id] != -float("inf"), "Clean codon was incorrectly masked!"


def test_constrained_beam_search_engine_optimize():
    engine = ConstrainedBeamSearchEngine()
    result = engine.optimize_cds(
        "MAKW", host="nbenthamiana", gc_band="40-47", terminal_stop_policy="append"
    )

    assert result["generation_engine"] == "lm"
    assert result["inference_mode"] == "deterministic_scaffold"
    assert result["trained_model_used"] is False
    assert result["optimized_sequence"].startswith("ATG")
    assert result["optimized_sequence"].endswith(("TAA", "TAG", "TGA"))
    assert translate_dna(result["optimized_sequence"]) == "MAKW*"
    report = EvaluationResult.model_validate(result["evaluation_report"])
    assert report.sequence_integrity.aa_identity == 1.0
    assert report.sequence_integrity.internal_stop_count == 0
    checks = {check.check_name: check.result.value for check in report.checks}
    assert checks["type_iis_bsai"] == "pass"
    assert checks["type_iis_bsmbi"] == "pass"
    assert checks["type_iis_bpii"] == "pass"
    # This deterministic short fixture falls below the requested GC band.
    # Translation validity must not turn its failed GC gate into an overall pass.
    assert report.metrics.gc_percent < 40
    assert checks["gc_content"] == "fail"
    assert result["validator_passed"] is report.passed is False
