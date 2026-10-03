"""Unit tests for Job 350: Experimental Codon UI & Watermark Pairs."""

import pytest
from factorforge.codon import (
    CodonDistributionPolicy,
    DOUG_BALANCED_PRESET,
    ILLUSTRATIVE_SKEWED_PRESET,
    DeterministicAllocator,
    TwinCandidateGenerator,
)
from factorforge.codon.qa import InSilicoQA, translate


LC_CDS_REF = (
    "ATGGCAAAAACAAATTTATTTTTATTTTTAATATTTTCATTATTACTGTCTCTGTCCTCTGCTGACATTCAAATGACTCAAAGCCCTTCTTCTCTTTCT"
    "GCTTCTGTTGGAGACAGAGTCACCATCACTTGCAGAGCTTCTCAAGGAATCAGAAATTACCTTGCTTGGTATCAACAGAAGCCTGGAAAAGCTCCCAAA"
    "CTTTTGATTTATGCTGCTTCCACTCTTCAATCTGGAGTTCCATCCAGATTTTCTGGATCTGGATCTGGAACTGATTTCACTTTGACCATTTCTTCTTTG"
    "CAACCTGAAGATGTTGCCACTTACTATTGTCAAAGATACAACAGAGCACCATACACTTTTGGACAGGGAACCAAAGTTGAAATCAAGAGAACTGTTGCT"
    "GCTCCTTCTGTTTTCATTTTTCCTCCTTCTGATGAACAATTGAAATCTGGAACTGCCTCTGTTGTTTGTCTTTTGAACAATTTTTATCCCAGAGAAGCC"
    "AAAGTTCAATGGAAGGTTGACAATGCTCTTCAATCTGGAAACAGCCAGGAATCTGTCACTGAACAAGATTCCAAAGATTCCACTTATTCTCTCTCTTCC"
    "ACTTTGACACTTTCCAAGGCTGATTATGAGAAACACAAAGTTTATGCATGTGAAGTCACTCATCAAGGACTTTCTTCTCCTGTCACCAAATCATTCAAC"
    "AGAGGAGAATGT"
)

HC_CDS_REF = (
    "ATGGCAAAAACAAATTTATTTTTATTTTTAATATTTTCATTATTACTGTCTCTGTCCTCTGCTGAAGTGCAATTGGTTGAGTCTGGAGGAGGACTGGTT"
    "CAACCTGGAAGATCTTTGAGACTTTCTTGTGCTGCTTCTGGATTCACTTTTGATGATTATGCCATGCATTGGGTCAGGCAGGCTCCTGGAAAGGGACTT"
    "GAATGGGTTTCTGCCATCACCTGGAATTCAGGGCACATTGATTATGCTGACTCTGTTGAAGGAAGATTCACCATTTCCAGAGACAATGCCAAAAATTCT"
    "CTTTATCTTCAGATGAATTCTTTGAGAGCTGAGGACACTGCTGTTTATTATTGTGCCAAAGTTTCTTATCTTTCCACTGCCAGCTCTTTGGATTATTGG"
    "GGACAAGGAACTCTTGTCACTGTTTCTTCTGCTTCCACAAAAGGGCCTTCTGTTTTTCCTCTTGCTCCTTCATCCAAATCCACTTCTGGTGGAACTGCT"
    "GCTCTTGGATGTCTTGTCAAGGATTATTTTCCTGAACCTGTCACTGTGTCTTGGAATTCTGGAGCTTTGACTTCTGGAGTTCACACATTTCCTGCAGTG"
    "CTTCAGTCCTCTGGCCTTTACAGTCTTTCTTCTGTTGTCACTGTTCCATCTTCTTCTTTGGGAACTCAAACTTACATTTGCAATGTCAATCACAAACCA"
    "TCCAACACCAAAGTTGACAAAAAAGTTGAACCCAAGTCTTGTGACAAGACTCACACTTGTCCTCCTTGTCCTGCACCTGAACTTCTTGGAGGACCTTCT"
    "GTTTTTCTTTTCCCTCCCAAGCCAAAAGACACTTTGATGATTTCCAGAACTCCTGAAGTCACTTGTGTTGTTGTTGATGTTTCTCATGAAGATCCTGAA"
    "GTCAAATTCAACTGGTATGTTGATGGAGTGGAAGTTCACAATGCCAAAACCAAACCCAGAGAAGAACAATACAATTCCACATACAGGGTTGTTTCTGTT"
    "TTGACTGTTCTTCATCAGGATTGGTTGAATGGAAAAGAGTACAAATGCAAAGTTTCCAACAAAGCTCTTCCTGCTCCCATTGAAAAAACCATTTCCAAA"
    "GCCAAAGGACAACCCAGAGAGCCTCAAGTCTACACTCTTCCTCCTTCCAGAGATGAATTGACCAAGAATCAGGTTTCTTTGACCTGTCTTGTCAAAGGA"
    "TTTTATCCTTCTGACATTGCTGTTGAGTGGGAATCCAATGGACAACCTGAAAACAATTACAAAACCACTCCCCCTGTTCTTGACAGTGATGGATCTTTT"
    "TTTCTTTATTCCAAATTGACTGTTGACAAATCCAGATGGCAACAAGGCAATGTTTTTTCTTGCAGTGTCATGCATGAGGCTCTTCACAATCATTACACC"
    "CAAAAATCTCTTTCTTTGTCTCCAGGAAAA"
)


class TestJob350CodonPolicyAndAllocator:
    def test_policy_validation_success(self):
        policy = CodonDistributionPolicy(
            policy_id="test_policy_valid",
            distributions={"R": {"AGA": 0.5, "AGG": 0.5}},
        )
        assert policy.policy_id == "test_policy_valid"
        assert len(policy.canonical_sha256) == 64

    def test_policy_validation_rejects_invalid_sum(self):
        with pytest.raises(ValueError, match="must sum to 1.0"):
            CodonDistributionPolicy(
                policy_id="bad_sum",
                distributions={"R": {"AGA": 0.6, "AGG": 0.6}},
            )

    def test_policy_validation_rejects_wrong_amino_acid_codon(self):
        with pytest.raises(ValueError, match="does not code for amino acid"):
            CodonDistributionPolicy(
                policy_id="bad_codon",
                distributions={"R": {"AGA": 0.5, "TTT": 0.5}},  # TTT is Phe, not Arg
            )

    def test_largest_remainder_allocation_odd_residue(self):
        # 9 Arg residues with 50:50 target -> exactly 5 AGA, 4 AGG
        res = DeterministicAllocator.allocate("R", 9, {"AGA": 0.5, "AGG": 0.5})
        assert res.total_residues == 9
        assert res.achieved_counts == {"AGA": 5, "AGG": 4}
        assert res.achieved_fractions["AGA"] == round(5 / 9, 6)
        assert res.achieved_fractions["AGG"] == round(4 / 9, 6)

    def test_largest_remainder_allocation_divisible_serine(self):
        # 30 Ser residues with 20% each -> exactly 6 of each
        targets = {"TCT": 0.2, "TCC": 0.2, "TCA": 0.2, "AGT": 0.2, "AGC": 0.2}
        res = DeterministicAllocator.allocate("S", 30, targets)
        assert res.total_residues == 30
        for c in targets:
            assert res.achieved_counts[c] == 6

    def test_largest_remainder_hc_serine(self):
        # 55 Ser residues in mature HC with 20% each -> exactly 11 of each
        targets = {"TCT": 0.2, "TCC": 0.2, "TCA": 0.2, "AGT": 0.2, "AGC": 0.2}
        res = DeterministicAllocator.allocate("S", 55, targets)
        assert res.total_residues == 55
        for c in targets:
            assert res.achieved_counts[c] == 11

    def test_skewed_allocation(self):
        # 9 Arg residues with 70:30 target -> 6 AGA, 3 AGG
        res = DeterministicAllocator.allocate("R", 9, {"AGA": 0.7, "AGG": 0.3})
        assert res.achieved_counts == {"AGA": 6, "AGG": 3}


class TestJob350TwinCandidateGeneration:
    def test_lc_twin_candidate_generation_invariants(self):
        gen = TwinCandidateGenerator(watermark_key="Eijex_FactorForge_Watermark_2026")
        prot_lc = translate(LC_CDS_REF)

        pair = gen.generate("LC", prot_lc, LC_CDS_REF, DOUG_BALANCED_PRESET)

        # Invariant 1: Translation matches original protein
        assert translate(pair.cds_w0).rstrip("*") == prot_lc.rstrip("*")
        assert translate(pair.cds_w1).rstrip("*") == prot_lc.rstrip("*")

        # Invariant 2: 5' Signal peptide is 100% frozen
        assert pair.cds_w0[:60] == LC_CDS_REF[:60]
        assert pair.cds_w1[:60] == LC_CDS_REF[:60]

        # Invariant 3: Quotas are 100% IDENTICAL between W0 and W1 (Zero Drift)
        assert pair.quotas_identical is True
        assert pair.arg_counts_w0 == pair.arg_counts_w1 == {"AGA": 5, "AGG": 4, "CGA": 0, "CGC": 0, "CGG": 0, "CGT": 0}
        assert pair.ser_counts_w0 == pair.ser_counts_w1
        assert sum(pair.ser_counts_w0.values()) == 30

        # Invariant 4: Type IIS restriction sites clean in both W0 and W1
        assert pair.w0_qa.type_iis_clean is True
        assert pair.w1_qa.type_iis_clean is True

        # Invariant 5: Watermark detected in W1, not detected in W0
        assert pair.w0_detection["is_detected"] is False
        assert pair.w1_detection["is_detected"] is True
        assert pair.w1_detection["score"] > pair.w0_detection["score"]

        # Invariant 6: Synonymous nucleotide delta exists
        assert pair.delta_nt > 50
        assert pair.gc_delta_percent > 8.0
        assert pair.w0_qa.cai is None
        assert pair.w0_qa.cai_status == "unavailable"
        assert pair.w0_qa.cai_reason == "host_codon_weight_reference_not_provided"

    def test_hc_twin_candidate_generation_invariants(self):
        gen = TwinCandidateGenerator(watermark_key="Eijex_FactorForge_Watermark_2026")
        prot_hc = translate(HC_CDS_REF)

        pair = gen.generate("HC", prot_hc, HC_CDS_REF, DOUG_BALANCED_PRESET)

        # Invariant 1: Translation
        assert translate(pair.cds_w0).rstrip("*") == prot_hc.rstrip("*")
        assert translate(pair.cds_w1).rstrip("*") == prot_hc.rstrip("*")

        # Invariant 2: SP frozen
        assert pair.cds_w0[:60] == HC_CDS_REF[:60]
        assert pair.cds_w1[:60] == HC_CDS_REF[:60]

        # Invariant 3: Zero Quota Drift
        assert pair.quotas_identical is True
        assert pair.arg_counts_w0 == pair.arg_counts_w1
        assert pair.ser_counts_w0 == pair.ser_counts_w1

        # Invariant 4: Type IIS clean
        assert pair.w0_qa.type_iis_clean is True
        assert pair.w1_qa.type_iis_clean is True

        # Invariant 5: Watermark detectability
        assert pair.w0_detection["is_detected"] is False
        assert pair.w1_detection["is_detected"] is True
