# factorforge/tests/test_sllm_biological_pipeline.py
"""Comprehensive Test Suite for sLLM Biological Model Pipeline (Job 285A, 285B, 285C).

Validates:
1. CodonSequenceTokenizer invariants & synonymous masks.
2. GenomicCorpusBuilder QC filters & length alignment.
3. Statistical baselines (Baseline 0, 1, 2) evaluation.
4. CompactCodonTransformer architecture & synonymous masked loss.
5. OnnxSLMAdapter dynamic inference with quantized host model.
6. Constrained decoder multi-candidate slate generation under ResolvedDesignContract (100% pass rate).
"""

import pytest
import numpy as np

from factorforge.discovery.filter import HardConstraintConfig
from factorforge.engines.sllm.adapters import OnnxSLMAdapter
from factorforge.engines.sllm.compiler import DesignContractCompiler
from factorforge.engines.sllm.data.corpus_builder import GenomicCorpusBuilder
from factorforge.engines.sllm.data.tokenizer import (
    AA_TO_INDEX,
    CODON_TO_INDEX,
    CodonSequenceTokenizer,
)
from factorforge.engines.sllm.decoder import ConstrainedBeamDecoder, ConstrainedSampleDecoder
from factorforge.engines.sllm.model.baselines import (
    CodonBigramMarkovBaseline,
    CodonFrequencyBaseline,
)


def test_codon_tokenizer_invariants():
    tokenizer = CodonSequenceTokenizer()

    # Test codon encoding and decoding
    cds = "ATGGCATTTTAA"
    tokens = tokenizer.encode_codons(cds)
    assert len(tokens) == 4
    assert tokens[0] == CODON_TO_INDEX["ATG"]
    assert tokenizer.decode_codons(tokens) == cds

    # Test amino acid encoding and decoding
    aa = "MAF*"
    aa_tokens = tokenizer.encode_amino_acids(aa)
    assert len(aa_tokens) == 3  # stop codon stripped from AA sequence
    assert tokenizer.decode_amino_acids(aa_tokens) == "MAF"

    # Test synonymous mask
    ala_id = AA_TO_INDEX["A"]
    ala_mask = tokenizer.get_synonymous_mask(ala_id)
    assert ala_mask.shape == (64,)
    assert ala_mask.dtype == bool
    assert np.sum(ala_mask) == 4  # GCT, GCC, GCA, GCG


def test_genomic_corpus_builder_qc():
    builder = GenomicCorpusBuilder(min_aa_len=10, max_aa_len=50)

    # Valid record
    valid_dna = "ATG" + "GCT" * 15 + "TAA"
    is_valid, reason, coding_dna, term_stop = builder.validate_transcript_qc(valid_dna)
    assert is_valid is True
    assert reason is None
    assert term_stop == "TAA"
    assert coding_dna == "ATG" + "GCT" * 15

    # Non-ATG start
    bad_dna = "GTG" + "GCT" * 15 + "TAA"
    is_valid, reason, _, _ = builder.validate_transcript_qc(bad_dna)
    assert is_valid is False
    assert "start" in reason.lower()


def test_statistical_baselines():
    toy_records = [
        {"coding_cds": "ATGGCTTTTGAA", "protein_aa": "MAFE"},
        {"coding_cds": "ATGGCCGATGAG", "protein_aa": "MADE"},
    ]

    # Baseline 0 & 1
    freq_bl = CodonFrequencyBaseline()
    freq_bl.fit(toy_records)
    nll_m, ppl_m = freq_bl.evaluate_nll(toy_records, mode="marginal")
    nll_c, ppl_c = freq_bl.evaluate_nll(toy_records, mode="conditional")
    assert nll_m > 0.0
    assert nll_c > 0.0
    assert ppl_m >= 1.0

    # Baseline 2
    bigram_bl = CodonBigramMarkovBaseline()
    bigram_bl.fit(toy_records)
    nll_b2, ppl_b2 = bigram_bl.evaluate_nll(toy_records)
    assert nll_b2 > 0.0


def test_compact_codon_transformer_forward():
    torch = pytest.importorskip("torch")
    from factorforge.engines.sllm.model.arch import (
        CompactCodonTransformer,
        compute_synonymous_masked_loss,
    )

    model = CompactCodonTransformer(
        d_model=64,
        nhead=2,
        num_layers=2,
        dim_feedforward=128,
        max_len=128,
    )

    B, L = 2, 10
    input_codons = torch.randint(0, 64, (B, L))
    aa_ids = torch.randint(0, 20, (B, L))

    logits = model(input_codons, aa_ids)
    assert logits.shape == (B, L, 64)

    # Synonymous masked loss test
    syn_masks = torch.ones((B, L, 64), dtype=torch.bool)
    pad_mask = torch.ones((B, L), dtype=torch.bool)
    target_codons = torch.randint(0, 64, (B, L))

    loss = compute_synonymous_masked_loss(logits, target_codons, syn_masks, pad_mask)
    assert loss.item() > 0.0
    assert not torch.isnan(loss)


def test_onnx_model_adapter_inference():
    pytest.importorskip("onnxruntime")
    adapter = OnnxSLMAdapter()
    assert "onnx" in adapter.model_name

    logits = adapter.predict_logits(
        prefix_codons=["ATG", "GCT"],
        target_aa="K",
        position=2,
    )
    assert logits.scores.shape == (64,)
    assert not np.any(np.isnan(logits.scores))


def test_biological_slate_constrained_generation():
    pytest.importorskip("onnxruntime")
    contract = DesignContractCompiler.resolve_contract(
        host="nbenthamiana",
        config=HardConstraintConfig(
            forbidden_type_iis_enzymes={"BsaI", "BsmBI"},
            forbidden_other_motifs=["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"],
            homopolymer_max_run=5,
            global_gc_min=30.0,
            global_gc_max=60.0,
        ),
        expected_stop="TAA",
    )

    adapter = OnnxSLMAdapter()
    beam_decoder = ConstrainedBeamDecoder(contract=contract, model_adapter=adapter, beam_width=3)
    sample_decoder = ConstrainedSampleDecoder(
        contract=contract, model_adapter=adapter, temperature=0.7, seed=42
    )

    target_aa = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"

    cand_beam = beam_decoder.decode(target_aa, stop_codon="TAA")
    assert cand_beam["status"] in ("SUCCESS", "DP_RESCUED")

    cand_sample = sample_decoder.decode(target_aa, stop_codon="TAA")
    assert cand_sample["status"] in ("SUCCESS", "DP_RESCUED")

    # Verify independent validator passes 100%
    for cand in [cand_beam, cand_sample]:
        val = contract.final_validator.evaluate(
            sequence_dna=cand["sequence"],
            expected_protein_aa=target_aa,
            expected_stop="TAA",
        )
        assert val.is_feasible is True, f"Violations: {val.violations}"
