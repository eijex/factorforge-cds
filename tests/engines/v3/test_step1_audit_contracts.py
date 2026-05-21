"""Focused STEP 1 audit tests for metrics, mappings, and v2/v3 contracts."""

from __future__ import annotations

import pytest

from factorforge.engines.v2.optimizer import RuleBasedOptimizer
from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile, ReverseTranslator
from factorforge.engines.v3.metrics import (
    CodonUsageTable,
    compute_cai,
    compute_gc,
    load_codon_usage_table,
)
from factorforge.engines.v3.tokenizer import CODON_TOKENS, SPECIAL_TOKENS, CodonTokenizer


STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate_dna(dna_sequence: str, codon_to_aa: dict[str, str]) -> str:
    codons = [dna_sequence[index : index + 3] for index in range(0, len(dna_sequence), 3)]
    return "".join(codon_to_aa[codon] for codon in codons)


def test_gc_calculation_known_sequences() -> None:
    assert compute_gc("") == 0.0
    assert compute_gc("ATGC") == 50.0
    assert compute_gc("GGCC") == 100.0
    assert compute_gc("AATT") == 0.0


def test_cai_calculation_simple_known_table() -> None:
    table = CodonUsageTable(
        codon_to_aa={"ATG": "M", "GCC": "A", "GCT": "A"},
        codon_weights={"ATG": 1.0, "GCC": 1.0, "GCT": 0.25},
        best_codon_for_aa={"M": "ATG", "A": "GCC"},
    )

    assert compute_cai("ATGGCC", table) == pytest.approx(1.0)
    assert compute_cai("ATGGCT", table) == pytest.approx(0.5)


def test_dna_translation_known_sequence() -> None:
    table = load_codon_usage_table()
    assert translate_dna("ATGGCTTTC", table.codon_to_aa) == "MAF"


def test_codon_to_amino_acid_mapping_contains_standard_codons() -> None:
    table = load_codon_usage_table()
    assert table.codon_to_aa["ATG"] == "M"
    assert table.codon_to_aa["TGG"] == "W"
    assert {table.codon_to_aa[codon] for codon in STOP_CODONS} == {"*"}
    assert set(CODON_TOKENS) == set(table.codon_to_aa)


def test_reverse_translation_preserves_amino_acid_sequence() -> None:
    protein = "ACDEFGHIKLMNPQRSTVWY"
    translator = ReverseTranslator()
    table = load_codon_usage_table()

    dna = translator.reverse_translate(protein, profile=OptimizationProfile.HIGH_CAI)

    assert len(dna) == len(protein) * 3
    assert translate_dna(dna, table.codon_to_aa) == protein


def test_v2_output_has_no_internal_stop_codons() -> None:
    protein = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"
    result = RuleBasedOptimizer().optimize(protein, profile="high_cai")
    codons = [result.sequence[index : index + 3] for index in range(0, len(result.sequence), 3)]

    assert not any(codon in STOP_CODONS for codon in codons[:-1])


def test_every_standard_amino_acid_maps_only_to_valid_synonymous_codons() -> None:
    translator = ReverseTranslator()
    table = load_codon_usage_table()

    for aa, codons in translator.aa_to_codons.items():
        if aa == "*":
            continue
        assert codons, f"{aa} has no codons"
        for codon, _frequency in codons:
            assert codon in table.codon_to_aa
            assert table.codon_to_aa[codon] == aa
            assert codon not in STOP_CODONS


def test_tokenizer_special_tokens_are_not_codon_tokens_and_decode_skips_them() -> None:
    tokenizer = CodonTokenizer.default()
    assert set(SPECIAL_TOKENS).isdisjoint(CODON_TOKENS)

    encoded = [
        tokenizer.bos_token_id,
        tokenizer.token_to_id["ATG"],
        tokenizer.pad_token_id,
        tokenizer.token_to_id["GCC"],
        tokenizer.eos_token_id,
    ]

    assert tokenizer.decode(encoded, skip_special_tokens=True) == "ATGGCC"

