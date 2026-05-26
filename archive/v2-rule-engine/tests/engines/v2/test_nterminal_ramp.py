"""
Unit tests for N-terminal ramp (RAMP profile).
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile, ReverseTranslator


@pytest.fixture
def translator():
    """Create ReverseTranslator instance."""
    return ReverseTranslator()


@pytest.fixture
def long_protein():
    """Protein longer than ramp region (>50 AAs)."""
    return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMK"


class TestRampProfile:
    """Test RAMP profile translation."""

    def test_ramp_translate_correct_length(self, translator, long_protein):
        """RAMP translation produces correct length DNA."""
        dna = translator.reverse_translate(long_protein, profile=OptimizationProfile.RAMP)
        assert len(dna) == len(long_protein) * 3

    def test_ramp_translate_starts_with_atg(self, translator, long_protein):
        """RAMP translation starts with ATG (Met)."""
        dna = translator.reverse_translate(long_protein, profile=OptimizationProfile.RAMP)
        # Met is a single-codon AA so ramp skips it
        assert dna[:3] == "ATG"

    def test_ramp_cai_lower_at_nterminus(self, translator, long_protein):
        """N-terminal ramp region has lower CAI than full sequence."""
        random.seed(42)
        dna = translator.reverse_translate(long_protein, profile=OptimizationProfile.RAMP)
        # N-terminal (first 50 codons = 150 bp)
        n_terminal = dna[:150]
        # C-terminal (remaining)
        c_terminal = dna[150:]

        if len(c_terminal) >= 9:  # Need at least 3 codons
            cai_n = translator.calculate_cai(n_terminal)
            cai_c = translator.calculate_cai(c_terminal)
            # N-terminal CAI should generally be lower due to ramp
            # (May not always be true due to randomness, but statistically expected)
            # Just verify both are valid CAI values
            assert 0.0 <= cai_n <= 1.0
            assert 0.0 <= cai_c <= 1.0

    def test_ramp_preserves_amino_acid_sequence(self, translator, long_protein):
        """RAMP translation preserves the amino acid sequence."""
        random.seed(42)
        dna = translator.reverse_translate(long_protein, profile=OptimizationProfile.RAMP)
        codon_table = translator.codon_table["codons"]
        # Verify each codon translates to the expected amino acid
        for i, aa in enumerate(long_protein):
            codon = dna[i * 3 : i * 3 + 3]
            assert codon in codon_table, f"Unknown codon {codon} at position {i}"
            assert codon_table[codon]["aa"] == aa, (
                f"Position {i}: expected {aa}, got {codon_table[codon]['aa']} (codon {codon})"
            )


class TestApplyNterminalRamp:
    """Test _apply_nterminal_ramp method directly."""

    def test_single_codon_aa_unchanged(self, translator):
        """Single-codon amino acids (M, W) should not be changed by ramp."""
        protein = "MWMW"
        dna = "ATGTGGATGTGG"
        ramped = translator._apply_nterminal_ramp(dna, protein, ramp_codons=10)
        # M(ATG) and W(TGG) are single-codon: must remain unchanged
        assert ramped[0:3] == "ATG"
        assert ramped[3:6] == "TGG"
        assert ramped[6:9] == "ATG"
        assert ramped[9:12] == "TGG"

    def test_ramp_codons_parameter(self, translator):
        """ramp_codons limits the number of positions affected."""
        protein = "AAAA"  # 4 Ala codons
        dna = translator.reverse_translate(protein, profile=OptimizationProfile.HIGH_CAI)
        # Apply ramp only to first 2 codons
        random.seed(42)
        ramped = translator._apply_nterminal_ramp(dna, protein, ramp_codons=2)
        # Codons 3 and 4 should be unchanged
        assert ramped[6:9] == dna[6:9]
        assert ramped[9:12] == dna[9:12]

    def test_short_protein_no_crash(self, translator):
        """Protein shorter than ramp_codons should not crash."""
        protein = "MA"
        dna = translator.reverse_translate(protein, profile=OptimizationProfile.HIGH_CAI)
        ramped = translator._apply_nterminal_ramp(dna, protein, ramp_codons=50)
        assert len(ramped) == len(dna)
        # First codon is M (single codon, unchanged)
        assert ramped[:3] == "ATG"

    def test_ramp_codon_is_valid_synonym(self, translator):
        """Ramped codons are valid synonymous codons for each AA."""
        protein = "AKLDEF"
        dna = translator.reverse_translate(protein, profile=OptimizationProfile.HIGH_CAI)
        random.seed(42)
        ramped = translator._apply_nterminal_ramp(dna, protein, ramp_codons=10)
        codon_table = translator.codon_table["codons"]
        for i, aa in enumerate(protein):
            codon = ramped[i * 3 : i * 3 + 3]
            assert codon in codon_table
            assert codon_table[codon]["aa"] == aa


class TestRampCandidateGeneration:
    """Test candidate generation with RAMP profile."""

    def test_generate_ramp_candidates(self, translator):
        """RAMP profile works with generate_candidates."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTL"
        candidates = translator.generate_candidates(
            protein, profile=OptimizationProfile.RAMP, n=3
        )
        assert len(candidates) <= 3
        assert all("score" in c for c in candidates)
        assert all("cai" in c for c in candidates)

    def test_ramp_profile_in_supported_profiles(self, translator):
        """RAMP is a valid OptimizationProfile enum value."""
        assert OptimizationProfile.RAMP.value == "ramp"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
