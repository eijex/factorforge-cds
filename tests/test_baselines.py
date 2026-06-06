"""Baseline CDS generator tests — benchmark foundation."""
from factorforge.analysis.metrics import translate_dna
from benchmarks.baselines.random_synonymous import random_synonymous_cds


def test_random_synonymous_preserves_protein():
    protein = "MKTAY"
    cds = random_synonymous_cds(protein, seed=320)
    assert translate_dna(cds).rstrip("*") == protein


def test_random_synonymous_deterministic_with_seed():
    assert random_synonymous_cds("MKTAY", seed=320) == random_synonymous_cds("MKTAY", seed=320)


from benchmarks.baselines.most_frequent_codon import most_frequent_codon_cds


def test_most_frequent_preserves_protein():
    from factorforge.analysis.metrics import translate_dna
    assert translate_dna(most_frequent_codon_cds("MKTAY")).rstrip("*") == "MKTAY"


from benchmarks.baselines.greedy_cai import greedy_cai_cds
from factorforge.analysis.metrics import calculate_cai, load_codon_usage_table


def test_greedy_cai_high_cai():
    w = load_codon_usage_table().codon_weights
    cds = greedy_cai_cds("MKTAYIAK")
    # greedy should reach near-max CAI
    assert calculate_cai(cds, w) >= 0.9
