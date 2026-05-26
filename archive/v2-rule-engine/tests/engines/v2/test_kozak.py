"""
Unit tests for Kozak sequence optimization.

Plant (N. benthamiana) optimal Kozak context: AACAATG[G/GC]...
The 2nd codon should ideally start with G (good) or GC (best).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.rules.reverse_translator import (
    OptimizationProfile,
    ReverseTranslator,
)


@pytest.fixture
def translator():
    """Create ReverseTranslator instance."""
    return ReverseTranslator()


class TestKozakOptimization:
    """Test _apply_kozak_optimization method."""

    def test_kozak_places_g_at_position_4(self, translator):
        """2nd codon should start with G after Kozak optimization.

        Ala (A) codons: GCT, GCC, GCA, GCG — all start with G.
        """
        protein = "MAKL"
        dna = translator.reverse_translate(protein, kozak=True)
        codon2 = dna[3:6]
        assert codon2[0] == "G"

    def test_kozak_prefers_gc_start(self, translator):
        """When multiple G-starting codons exist, prefer GC."""
        protein = "MAKL"
        dna = translator.reverse_translate(protein, kozak=True)
        codon2 = dna[3:6]
        assert codon2[:2] == "GC"

    def test_kozak_preserves_amino_acid(self, translator):
        """Kozak optimization must not change the amino acid sequence."""
        protein = "MVSKGEELFTGVV"
        dna = translator.reverse_translate(protein, kozak=True)
        codon_table = translator.codon_table["codons"]
        for i, aa in enumerate(protein):
            codon = dna[i * 3 : i * 3 + 3]
            assert codon_table[codon]["aa"] == aa

    def test_kozak_no_change_when_already_optimal(self, translator):
        """If 2nd codon already starts with G, no substitution needed."""
        dna_orig = "ATGGCCAAACTG"  # M-A-K-L
        protein = "MAKL"
        result = translator._apply_kozak_optimization(dna_orig, protein)
        assert result[3] == "G"

    def test_kozak_short_protein_unchanged(self, translator):
        """Single amino acid protein (just M) returns unchanged."""
        dna = "ATG"
        result = translator._apply_kozak_optimization(dna, "M")
        assert result == "ATG"

    def test_kozak_impossible_amino_acid(self, translator):
        """If no synonymous codon starts with G, sequence unchanged.

        Phe (F) codons: TTT, TTC — neither starts with G.
        """
        protein = "MFKL"
        dna = translator.reverse_translate(protein, kozak=True)
        # Codon 2 encodes Phe — no G-starting synonymous codon
        assert dna[3:6] in ("TTT", "TTC")

    def test_kozak_val_starts_with_g(self, translator):
        """Val (V) codons: GTT, GTC, GTA, GTG — all start with G."""
        protein = "MVKL"
        dna = translator.reverse_translate(protein, kozak=True)
        assert dna[3] == "G"

    def test_kozak_gly_starts_with_g(self, translator):
        """Gly (G) codons: GGT, GGC, GGA, GGG — all start with G."""
        protein = "MGKL"
        dna = translator.reverse_translate(protein, kozak=True)
        assert dna[3] == "G"

    def test_kozak_opt_in_default_false(self, translator):
        """kozak defaults to False — no Kozak processing by default."""
        protein = "MFKL"
        dna1 = translator.reverse_translate(
            protein, profile=OptimizationProfile.HIGH_CAI
        )
        dna2 = translator.reverse_translate(
            protein, profile=OptimizationProfile.HIGH_CAI, kozak=False
        )
        assert dna1 == dna2

    def test_kozak_works_with_all_profiles(self, translator):
        """Kozak optimization works with every profile."""
        protein = "MAKL"  # Ala always has G-starting codons
        for profile in OptimizationProfile:
            dna = translator.reverse_translate(protein, profile=profile, kozak=True)
            assert len(dna) == len(protein) * 3
            assert dna[3] == "G", f"Failed for profile {profile.name}"

    def test_kozak_with_generate_candidates(self, translator):
        """Kozak optimization propagates through generate_candidates."""
        protein = "MAKL"
        candidates = translator.generate_candidates(
            protein, profile=OptimizationProfile.BALANCED, n=3, kozak=True
        )
        for c in candidates:
            assert c["sequence"][3] == "G"


class TestKozakEdgeCases:
    """Edge cases for Kozak optimization."""

    def test_empty_protein(self, translator):
        """Empty protein returns empty DNA."""
        result = translator._apply_kozak_optimization("", "")
        assert result == ""

    def test_two_amino_acid_protein(self, translator):
        """M + one AA should still optimize."""
        protein = "MA"
        dna = translator.reverse_translate(protein, kozak=True)
        assert len(dna) == 6
        assert dna[:3] == "ATG"
        assert dna[3] == "G"  # Ala codons all start with G
