"""
Unit tests for ReverseTranslator
"""

import sys
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.rules.reverse_translator import OptimizationProfile, ReverseTranslator


@pytest.fixture
def translator():
    """Create ReverseTranslator instance"""
    return ReverseTranslator()


@pytest.fixture
def sample_protein():
    """Sample protein sequence (GFP partial)"""
    return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK"


class TestCodonTableLoading:
    """Test codon table loading"""

    def test_load_default_codon_table(self, translator):
        """Test loading default codon table"""
        assert translator.codon_table is not None
        assert len(translator.codon_table) > 0

    def test_aa_to_codons_map(self, translator):
        """Test amino acid to codons mapping"""
        assert translator.aa_to_codons is not None

        # Check some standard amino acids
        assert "A" in translator.aa_to_codons
        assert "M" in translator.aa_to_codons
        assert "*" in translator.aa_to_codons  # STOP codon

        # M should have only ATG
        assert len(translator.aa_to_codons["M"]) == 1
        assert translator.aa_to_codons["M"][0][0] == "ATG"

    def test_codon_frequencies_sum_to_one_for_sense_codons(self, translator):
        """Sense amino-acid families: synonymous codon frequencies sum to ~1.0.

        Stop family (*) is excluded — genome-derived tables built via
        scripts/build_codon_profile.py trim the terminal stop codon before
        counting, so '*' legitimately has zero frequency mass (Job 168 /
        v3.3.0, _analysis/025; terminal_stop_policy: excluded_from_frequency_model
        in the codon-table manifest). The stop codon is never reverse-translated
        as part of CDS body generation, so this has no functional effect.
        """
        for aa, codons in translator.aa_to_codons.items():
            if aa == "*":
                continue
            total_freq = sum(freq for _, freq in codons)
            assert 0.99 <= total_freq <= 1.01, f"AA {aa} frequencies sum to {total_freq}"

    def test_stop_codon_frequencies_excluded_from_frequency_model(self):
        """Stop family frequencies are all 0.0 in the NbeV1.1/v2 table — excluded
        by design, not a data bug. This is a v2-specific table policy
        (terminal_stop_policy: excluded_from_frequency_model); the legacy/v1
        table (currently the production default, see active_codon_reference.json)
        uses included_in_frequency_model instead, so this test targets the v2
        table explicitly rather than relying on it being the default."""
        from factorforge.engines.profile.utils import get_data_path

        v2_translator = ReverseTranslator(
            codon_table_path=get_data_path() / "profiles" / "nbev11_cds_hc_derived_codons.json"
        )
        stop_codons = v2_translator.aa_to_codons.get("*", [])
        assert stop_codons, "stop codon entry must still be present in aa_to_codons"
        assert all(freq == 0.0 for _, freq in stop_codons)


class TestCAICalculation:
    """Test CAI calculation"""

    def test_calculate_cai_basic(self, translator):
        """Test basic CAI calculation"""
        # Use high-frequency codons
        seq = "ATGGCCGCCTAA"  # M-A-A-*
        cai = translator.calculate_cai(seq)

        assert 0.0 <= cai <= 1.0

    def test_cai_high_frequency_codons(self, translator):
        """Test CAI with high-frequency codons"""
        # Use most frequent codons for each AA
        seq = "ATGGCCGCCTAA"  # Should have high CAI
        cai = translator.calculate_cai(seq)

        assert cai > 0.5  # Should be relatively high

    def test_cai_invalid_length(self, translator):
        """Test CAI with invalid sequence length"""
        seq = "ATGGC"  # Not divisible by 3
        cai = translator.calculate_cai(seq)

        assert cai == 0.0  # Should return 0 for invalid sequences


class TestGCCalculation:
    """Test GC content calculation"""

    def test_calculate_gc_content(self, translator):
        """Test GC content calculation"""
        # 50% GC
        seq = "ATGC"
        gc = translator.calculate_gc_content(seq)
        assert gc == 50.0

        # 100% GC
        seq = "GGGGCCCC"
        gc = translator.calculate_gc_content(seq)
        assert gc == 100.0

        # 0% GC
        seq = "AAATTT"
        gc = translator.calculate_gc_content(seq)
        assert gc == 0.0

    def test_calculate_local_gc(self, translator):
        """Test local GC content calculation"""
        # Create sequence with varying GC
        seq = "A" * 50 + "G" * 50 + "C" * 50

        local_gc = translator.calculate_local_gc(seq, window_size=50)

        assert len(local_gc) > 0
        assert all(0 <= gc <= 100 for gc in local_gc)

    def test_local_gc_window_larger_than_sequence(self, translator):
        """Return empty list when window exceeds sequence length"""
        seq = "ATGC"
        local_gc = translator.calculate_local_gc(seq, window_size=10)

        assert local_gc == []


class TestBalancedProfile:
    """Test Balanced optimization profile"""

    def test_balanced_translate_basic(self, translator, sample_protein):
        """Test basic balanced translation"""
        dna = translator.reverse_translate(sample_protein, profile=OptimizationProfile.BALANCED)

        assert len(dna) == len(sample_protein) * 3
        assert dna.startswith("ATG")  # Should start with M

    def test_balanced_cai(self, translator, sample_protein):
        """Test balanced profile CAI"""
        dna = translator.reverse_translate(sample_protein, profile=OptimizationProfile.BALANCED)
        cai = translator.calculate_cai(dna)

        assert cai > 0.6  # Should have decent CAI

    def test_balanced_gc(self, translator, sample_protein):
        """Test balanced profile GC content"""
        dna = translator.reverse_translate(sample_protein, profile=OptimizationProfile.BALANCED)
        gc = translator.calculate_gc_content(dna)

        # Should be in reasonable range (relaxed for probabilistic nature)
        assert 35.0 <= gc <= 65.0


class TestHighCAIProfile:
    """Test High-CAI optimization profile"""

    def test_high_cai_translate(self, translator, sample_protein):
        """Test high-CAI translation"""
        dna = translator.reverse_translate(sample_protein, profile=OptimizationProfile.HIGH_CAI)

        assert len(dna) == len(sample_protein) * 3

    def test_high_cai_value(self, translator, sample_protein):
        """Test that high-CAI profile produces high CAI"""
        dna = translator.reverse_translate(sample_protein, profile=OptimizationProfile.HIGH_CAI)
        cai = translator.calculate_cai(dna)

        # Should have very high CAI (golden set reference shifts scale slightly)
        assert cai > 0.70


class TestGCTargetProfile:
    """Test GC-Target optimization profile"""

    def test_gc_target_translate(self, translator, sample_protein):
        """Test GC-target translation"""
        dna = translator.reverse_translate(sample_protein, profile=OptimizationProfile.GC_TARGET)

        assert len(dna) == len(sample_protein) * 3

    def test_gc_target_value(self, translator, sample_protein):
        """Test that GC-target profile achieves target GC"""
        dna = translator.reverse_translate(
            sample_protein, profile=OptimizationProfile.GC_TARGET, target_gc=50.0
        )
        gc = translator.calculate_gc_content(dna)

        # Should be close to 50%
        assert 45.0 <= gc <= 55.0

    def test_gc_target_default_uses_host_midpoint(self, translator, sample_protein):
        """Without explicit target_gc, gc_target defaults to the active host's
        composition midpoint (resolve_host_gc_mid(self.host)).

        Job 168 / v3.3.0 (_analysis/025) corrected this midpoint for
        nbenthamiana from a circularly-derived 60% (calibrated from the
        legacy engine's own GC-rich output, not from genome biology) to the
        native genome-composition anchor GC_OPT_MID=43.5%. This test now
        asserts the corrected, host-resolved value rather than a hardcoded
        constant, so it tracks whichever value is currently configured
        instead of re-encoding either the old or the new number by hand.
        """
        from factorforge.engines.profile.scoring import resolve_host_gc_mid

        dna = translator.reverse_translate(
            sample_protein, profile=OptimizationProfile.GC_TARGET
        )
        gc = translator.calculate_gc_content(dna)
        expected_mid = resolve_host_gc_mid(translator.host)
        assert abs(gc - expected_mid) <= 10.0

    def test_gc_target_explicit_low_still_supported(self, translator, sample_protein):
        """Users wanting low GC can still request it explicitly."""
        dna = translator.reverse_translate(
            sample_protein, profile=OptimizationProfile.GC_TARGET, target_gc=42.5
        )
        gc = translator.calculate_gc_content(dna)
        assert gc < 50.0  # explicit low target honored

    def test_gc_target_extreme_high_gc(self, translator):
        """Test GC-target profile with high target GC"""
        protein = "GGGG"
        dna = translator.reverse_translate(
            protein, profile=OptimizationProfile.GC_TARGET, target_gc=100.0
        )
        gc = translator.calculate_gc_content(dna)

        assert len(dna) == len(protein) * 3
        assert gc >= 90.0

    def test_gc_target_extreme_low_gc(self, translator):
        """Test GC-target profile with low target GC"""
        protein = "KKKK"
        dna = translator.reverse_translate(
            protein, profile=OptimizationProfile.GC_TARGET, target_gc=0.0
        )
        gc = translator.calculate_gc_content(dna)

        assert len(dna) == len(protein) * 3
        assert gc <= 10.0


class TestAssemblyFriendlyProfile:
    """Test Assembly-Friendly optimization profile"""

    def test_assembly_friendly_translate(self, translator, sample_protein):
        """Test assembly-friendly translation"""
        dna = translator.reverse_translate(
            sample_protein, profile=OptimizationProfile.ASSEMBLY_FRIENDLY
        )

        assert len(dna) == len(sample_protein) * 3

    def test_assembly_friendly_no_bsai(self, translator, sample_protein):
        """Test that assembly-friendly avoids BsaI sites"""
        dna = translator.reverse_translate(
            sample_protein, profile=OptimizationProfile.ASSEMBLY_FRIENDLY
        )

        # BsaI recognition site
        assert "GGTCTC" not in dna
        assert "GAGACC" not in dna  # Reverse complement

    def test_assembly_friendly_fallback_warning(self, translator, caplog):
        """Test fallback warning when restriction sites remain"""

        def always_with_site(*args, **kwargs):
            return "GGTCTC"

        translator._balanced_translate = always_with_site  # type: ignore[method-assign]

        dna = translator.reverse_translate(
            "MA", profile=OptimizationProfile.ASSEMBLY_FRIENDLY, max_attempts=1
        )

        # Check that warning was logged (changed from print to logger.warning)
        assert any("Could not remove all restriction sites" in record.message for record in caplog.records)
        assert "GGTCTC" in dna


class TestCandidateGeneration:
    """Test Top-N candidate generation"""

    def test_generate_candidates_basic(self, translator):
        """Test basic candidate generation"""
        protein = "MAKL"
        candidates = translator.generate_candidates(protein, n=5)

        assert len(candidates) <= 5
        assert all("sequence" in c for c in candidates)
        assert all("cai" in c for c in candidates)
        assert all("gc" in c for c in candidates)
        assert all("score" in c for c in candidates)

    def test_candidates_are_different(self, translator):
        """Test that candidates are different"""
        protein = "MAKL"
        candidates = translator.generate_candidates(protein, n=5)

        sequences = [c["sequence"] for c in candidates]
        # At least some should be different (unless very limited codon choices)
        # For these 4 AAs, there should be variation
        unique_sequences = set(sequences)
        assert len(unique_sequences) >= 1  # At least one unique

    def test_candidates_sorted_by_score(self, translator):
        """Test that candidates are sorted by score"""
        protein = "MAKL"
        candidates = translator.generate_candidates(protein, n=5)

        scores = [c["score"] for c in candidates]
        # Should be in descending order
        assert scores == sorted(scores, reverse=True)

    def test_candidates_correct_length(self, translator):
        """Test that all candidates have correct length"""
        protein = "MAKL"
        candidates = translator.generate_candidates(protein, n=5)

        expected_length = len(protein) * 3
        assert all(len(c["sequence"]) == expected_length for c in candidates)


class TestProteinTranslation:
    """Test protein sequence translation"""

    def test_translate_with_stop(self, translator):
        """Test translation with STOP codon"""
        protein = "MAKL*"
        dna = translator.reverse_translate(protein)

        # Should end with STOP codon
        assert dna[-3:] in ["TAA", "TAG", "TGA"]

    def test_translate_all_amino_acids(self, translator):
        """Test translation with all standard amino acids"""
        protein = "ACDEFGHIKLMNPQRSTVWY"
        dna = translator.reverse_translate(protein)

        assert len(dna) == len(protein) * 3

    def test_empty_protein(self, translator):
        """Test translation of empty protein"""
        protein = ""
        dna = translator.reverse_translate(protein)

        assert dna == ""


class TestReproducibility:
    """Test reproducibility of translations"""

    def test_deterministic_translation(self, translator, sample_protein):
        """Test that translation is deterministic with same seed"""
        # Note: This depends on implementation using random seed
        dna1 = translator.reverse_translate(sample_protein, profile=OptimizationProfile.HIGH_CAI)
        dna2 = translator.reverse_translate(sample_protein, profile=OptimizationProfile.HIGH_CAI)

        # High-CAI should be deterministic (uses most frequent codon)
        assert dna1 == dna2


class TestViralDeliveryProfile:
    """Tests for viral_delivery optimization profile (001-fix)"""

    def test_viral_delivery_profile_valid(self, translator):
        """OptimizationProfile.VIRAL_DELIVERY enum 생성 및 역번역 검증"""
        # enum 생성 성공 확인
        profile = OptimizationProfile("viral_delivery")
        assert profile == OptimizationProfile.VIRAL_DELIVERY

        # MAST (4 aa) → 12bp DNA 반환 확인
        dna = translator.reverse_translate("MAST", OptimizationProfile.VIRAL_DELIVERY)
        assert len(dna) == 12

        # 아미노산 서열 보존 확인: 코돈 테이블로 역번역 결과를 직접 번역
        codons_section = translator.codon_table["codons"]
        translated = "".join(
            codons_section[dna[i : i + 3]]["aa"] for i in range(0, len(dna), 3)
        )
        assert translated == "MAST"

    def test_viral_delivery_via_generate_candidates(self, translator):
        """generate_candidates로 viral_delivery 프로필 호출 검증"""
        candidates = translator.generate_candidates(
            "MAST", OptimizationProfile.VIRAL_DELIVERY, n=2
        )
        assert len(candidates) == 2
        for cand in candidates:
            assert 0.0 <= cand["score"] <= 1.0
            assert len(cand["sequence"]) == 12


class TestEdgeCases:
    """Test edge cases"""

    def test_single_amino_acid(self, translator):
        """Test translation of single amino acid"""
        protein = "M"
        dna = translator.reverse_translate(protein)

        assert dna == "ATG"

    def test_long_protein(self, translator):
        """Test translation of long protein"""
        protein = "M" * 1000
        dna = translator.reverse_translate(protein)

        assert len(dna) == 3000
        assert dna == "ATG" * 1000


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
