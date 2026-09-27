# factorforge/src/factorforge/engines/sllm/data/__init__.py
"""Data pipeline and tokenizers for FactorForge sLLM (Job 285A)."""

from factorforge.engines.sllm.data.corpus_builder import (
    CDSRecord,
    CorpusManifest,
    GenomicCorpusBuilder,
)
from factorforge.engines.sllm.data.tokenizer import (
    CodonSequenceTokenizer,
    CODON_VOCAB_SIZE,
    AA_VOCAB_SIZE,
    PAD_CODON_ID,
    BOS_CODON_ID,
    EOS_CODON_ID,
)

__all__ = [
    "CDSRecord",
    "CorpusManifest",
    "GenomicCorpusBuilder",
    "CodonSequenceTokenizer",
    "CODON_VOCAB_SIZE",
    "AA_VOCAB_SIZE",
    "PAD_CODON_ID",
    "BOS_CODON_ID",
    "EOS_CODON_ID",
]
