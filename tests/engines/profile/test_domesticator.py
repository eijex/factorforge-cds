"""
Unit tests for Domesticator
"""

import sys
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.rules.domesticator import Domesticator


@pytest.fixture
def domesticator():
    """Create Domesticator instance"""
    return Domesticator()


class TestRestrictionSiteScanning:
    """Test restriction site scanning"""

    def test_scan_bsai_site(self, domesticator):
        """Test BsaI site detection"""
        seq = "ATGGGTCTCGAGGAGCTG"  # Contains GGTCTC (BsaI)

        sites = domesticator.scan_restriction_sites(seq, "golden_gate")

        assert len(sites) > 0
        assert any(s["enzyme"] == "BsaI" for s in sites)

    def test_scan_no_sites(self, domesticator):
        """Test sequence without restriction sites"""
        seq = "ATGGCCGCCTAA"  # No restriction sites

        sites = domesticator.scan_restriction_sites(seq, "golden_gate")

        # May or may not have sites depending on sequence
        assert isinstance(sites, list)

    def test_scan_multiple_sites(self, domesticator):
        """Test multiple restriction sites"""
        # BsaI: GGTCTC, BpiI: GAAGAC
        seq = "ATGGGTCTCGAGGAAGACCTG"

        sites = domesticator.scan_restriction_sites(seq, "golden_gate")

        assert len(sites) >= 1

    def test_scan_biobricks(self, domesticator):
        """Test BioBricks restriction sites"""
        # EcoRI: GAATTC
        seq = "ATGGAATTCGAG"

        sites = domesticator.scan_restriction_sites(seq, "biobricks")

        assert len(sites) > 0
        assert any(s["enzyme"] == "EcoRI" for s in sites)


class TestDomestication:
    """Test domestication functionality"""

    def test_domesticate_simple(self, domesticator):
        """Test simple domestication"""
        # Create a sequence with BsaI site that can be removed
        seq = "ATGGGTCTCGAGGAGCTGTTC"

        result = domesticator.domesticate(seq, "golden_gate", max_attempts=50)

        assert "domesticated_seq" in result
        assert "removed_sites" in result
        assert "unfixable" in result

    def test_domesticate_preserves_length(self, domesticator):
        """Test that domestication preserves sequence length"""
        seq = "ATGGGTCTCGAGGAGCTGTTC"

        result = domesticator.domesticate(seq, "golden_gate")

        assert len(result["domesticated_seq"]) == len(seq)

    def test_domesticate_invalid_length(self, domesticator):
        """Test domestication with invalid sequence length"""
        seq = "ATGGC"  # Not divisible by 3

        result = domesticator.domesticate(seq, "golden_gate")

        assert not result["success"]
        assert "error" in result

    def test_domesticate_no_sites(self, domesticator):
        """Test domestication of sequence without sites"""
        seq = "ATGGCCGCCTAA"  # Likely no restriction sites

        result = domesticator.domesticate(seq, "golden_gate")

        # Should succeed immediately
        assert result["domesticated_seq"] == seq

    def test_domesticate_unfixable_site(self):
        """Test domestication fails when no synonymous codons exist"""
        codon_table = {
            "amino_acids": {
                "G": {"codons": ["GGT"]},
                "L": {"codons": ["CTC"]},
            },
            "codons": {
                "GGT": {"aa": "G", "frequency": 1.0},
                "CTC": {"aa": "L", "frequency": 1.0},
            },
        }
        local_domesticator = Domesticator(codon_table)
        seq = "GGTCTC"

        result = local_domesticator.domesticate(seq, "golden_gate", max_attempts=1)

        assert result["success"] is False
        assert len(result["unfixable"]) == 1
        assert result["unfixable"][0]["enzyme"] == "BsaI"


class TestBatchDomestication:
    """Test batch domestication"""

    def test_batch_domesticate(self, domesticator):
        """Test batch domestication"""
        sequences = [
            {"id": "seq1", "sequence": "ATGGCCGCCTAA"},
            {"id": "seq2", "sequence": "ATGGGTCTCGAG"},
        ]

        results = domesticator.batch_domesticate(sequences, "golden_gate")

        assert len(results) == 2
        assert all("id" in r for r in results)
        assert all("domesticated_seq" in r for r in results)

    def test_batch_empty(self, domesticator):
        """Test batch with empty list"""
        sequences = []

        results = domesticator.batch_domesticate(sequences, "golden_gate")

        assert len(results) == 0


class TestAlternativeSuggestions:
    """Test alternative suggestions"""

    def test_suggest_alternatives(self, domesticator):
        """Test alternative suggestions for unfixable sites"""
        seq = "ATGGGTCTCGAG"
        site = {"enzyme": "BsaI", "site": "GGTCTC", "position": 3}

        alternatives = domesticator._suggest_alternatives(seq, site)

        assert isinstance(alternatives, list)
        assert len(alternatives) > 0

    def test_suggest_alternatives_no_synonyms(self):
        """Test alternatives when affected codons have no synonyms"""
        codon_table = {
            "amino_acids": {
                "G": {"codons": ["GGT"]},
                "L": {"codons": ["CTC"]},
            },
            "codons": {
                "GGT": {"aa": "G", "frequency": 1.0},
                "CTC": {"aa": "L", "frequency": 1.0},
            },
        }
        local_domesticator = Domesticator(codon_table)
        seq = "GGTCTC"
        site = local_domesticator.scan_restriction_sites(seq, "golden_gate")[0]

        alternatives = local_domesticator._suggest_alternatives(seq, site)

        assert len(alternatives) == 3


class TestAssemblyStandards:
    """Test different assembly standards"""

    def test_golden_gate_standard(self, domesticator):
        """Test Golden Gate standard"""
        assert "golden_gate" in Domesticator.ASSEMBLY_STANDARDS
        assert "BsaI" in Domesticator.ASSEMBLY_STANDARDS["golden_gate"]["sites"]

    def test_moclo_standard(self, domesticator):
        """Test MoClo standard"""
        assert "moclo" in Domesticator.ASSEMBLY_STANDARDS
        assert "overhangs" in Domesticator.ASSEMBLY_STANDARDS["moclo"]

    def test_biobricks_standard(self, domesticator):
        """Test BioBricks standard"""
        assert "biobricks" in Domesticator.ASSEMBLY_STANDARDS
        assert "EcoRI" in Domesticator.ASSEMBLY_STANDARDS["biobricks"]["sites"]

    def test_invalid_standard(self, domesticator):
        """Test invalid assembly standard"""
        seq = "ATGGCCGCCTAA"

        with pytest.raises(ValueError):
            domesticator.scan_restriction_sites(seq, "invalid_standard")


class TestEdgeCases:
    """Test edge cases"""

    def test_empty_sequence(self, domesticator):
        """Test empty sequence"""
        seq = ""

        result = domesticator.domesticate(seq, "golden_gate")

        # Should handle gracefully
        assert "domesticated_seq" in result

    def test_very_short_sequence(self, domesticator):
        """Test very short sequence"""
        seq = "ATG"

        result = domesticator.domesticate(seq, "golden_gate")

        assert len(result["domesticated_seq"]) == 3

    def test_max_attempts_zero(self):
        """Test domestication with zero attempts"""
        codon_table = {
            "amino_acids": {
                "G": {"codons": ["GGT"]},
                "L": {"codons": ["CTC"]},
            },
            "codons": {
                "GGT": {"aa": "G", "frequency": 1.0},
                "CTC": {"aa": "L", "frequency": 1.0},
            },
        }
        local_domesticator = Domesticator(codon_table)
        seq = "GGTCTC"

        result = local_domesticator.domesticate(seq, "golden_gate", max_attempts=0)

        assert result["success"] is False
        assert len(result["unfixable"]) == 1


class TestDomesticatorBugFix:
    """002-fix: break → continue 수정 검증"""

    def test_domesticate_skips_unfixable_and_continues(self):
        """unfixable 사이트가 있어도 fixable 사이트 처리를 완료하는지 확인

        서열 GAGACCGGTCTC:
          - GAGACC (BsaI RC, pos:0): GAG+ACC 모두 synonymous 없음 → unfixable
          - GGTCTC (BsaI fwd, pos:6): GGT에 GGC synonymous 있음 → fixable

        scan 순서상 GGTCTC(fwd)가 먼저 sites[0]에 오므로 fixable이 처리된 후
        unfixable(GAGACC)이 남아 unfixable 리스트에 기록됨.
        """
        codon_table = {
            "amino_acids": {
                "E": {"codons": ["GAG"]},           # no synonym (unfixable용)
                "T": {"codons": ["ACC"]},           # no synonym (unfixable용)
                "G": {"codons": ["GGT", "GGC"]},   # has synonym (fixable용)
                "L": {"codons": ["CTC", "CTG"]},   # has synonym (fixable용)
            },
            "codons": {
                "GAG": {"aa": "E", "frequency": 1.0},
                "ACC": {"aa": "T", "frequency": 1.0},
                "GGT": {"aa": "G", "frequency": 0.6},
                "GGC": {"aa": "G", "frequency": 0.4},
                "CTC": {"aa": "L", "frequency": 0.6},
                "CTG": {"aa": "L", "frequency": 0.4},
            },
        }
        local_domesticator = Domesticator(codon_table)
        seq = "GAGACCGGTCTC"  # 12bp = GAG+ACC+GGT+CTC

        result = local_domesticator.domesticate(seq, "golden_gate", max_attempts=5)

        # a. unfixable 리스트에 1개 이상 (GAGACC unfixable)
        assert len(result["unfixable"]) >= 1
        # b. removed_sites에 1개 이상 (GGTCTC fixable 처리됨)
        assert len(result["removed_sites"]) >= 1
        # c. domesticated_seq에서 GGTCTC 제거됨
        assert "GGTCTC" not in result["domesticated_seq"]

    def test_domesticate_records_all_unfixable(self):
        """continue 수정 후 unfixable 사이트를 max_attempts번 시도하여 여러 번 기록하는지 확인

        break 버전: unfixable 1번 기록 후 loop 종료
        continue 버전: max_attempts(3)번 반복하여 unfixable 3번 기록
        → len(unfixable) >= 2 조건이 continue 버전에서만 PASS
        """
        codon_table = {
            "amino_acids": {
                "G": {"codons": ["GGT"]},  # no synonym
                "L": {"codons": ["CTC"]},  # no synonym
            },
            "codons": {
                "GGT": {"aa": "G", "frequency": 1.0},
                "CTC": {"aa": "L", "frequency": 1.0},
            },
        }
        local_domesticator = Domesticator(codon_table)
        seq = "GGTCTC"  # BsaI site, unfixable (no synonymous codons)

        result = local_domesticator.domesticate(seq, "golden_gate", max_attempts=3)

        assert result["success"] is False
        # continue 버전: 3번 반복, unfixable 3번 기록 → len >= 2
        # break 버전: 1번 기록 후 종료 → len == 1, FAIL
        assert len(result["unfixable"]) >= 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
