"""Baseline CDS generator tests — benchmark foundation."""
from factorforge.analysis.metrics import translate_dna
from benchmarks.baselines.random_synonymous import random_synonymous_cds


def test_random_synonymous_preserves_protein():
    protein = "MKTAY"
    cds = random_synonymous_cds(protein, seed=320)
    assert translate_dna(cds).rstrip("*") == protein


def test_random_synonymous_deterministic_with_seed():
    assert random_synonymous_cds("MKTAY", seed=320) == random_synonymous_cds("MKTAY", seed=320)
