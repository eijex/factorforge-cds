"""Integration tests for DiscoverySlateEngine end-to-end execution (Job 283A)."""

from factorforge.discovery.filter import HardConstraintFilter
from factorforge.discovery.slate import DiscoverySlateEngine
from factorforge.discovery.traits import TraitVectorExtractor


def test_hard_constraint_filter_validation():
    filt = HardConstraintFilter()
    target_aa = "MAW"

    # Valid DNA: M(ATG) A(GCT) W(TGG) *(TAA) -> ATGGCTTGATAA (or ATGGCTTGGTAA)
    # ATGGCTTGGTAA GC = 5/12 = 41.67%, len=12
    valid_dna = "ATGGCTTGGTAA"
    res = filt.evaluate(sequence_dna=valid_dna, expected_protein_aa=target_aa)
    assert res.is_feasible is True
    assert len(res.violations) == 0

    # Invalid: Contains BsaI (GGTCTC)
    invalid_bsai = "ATG" + "GGTCTC" + "TGG" + "TAA"  # M G L W *
    res_bsai = filt.evaluate(sequence_dna=invalid_bsai, expected_protein_aa="MGLW")
    assert res_bsai.is_feasible is False
    assert any("Forbidden restriction motifs" in v for v in res_bsai.violations)

    # Invalid: Contains NotI (GCGGCCGC)
    invalid_noti = "ATG" + "GCGGCCGC" + "TAA"
    res_noti = filt.evaluate(sequence_dna=invalid_noti, expected_protein_aa="MAA")
    assert res_noti.is_feasible is False


def test_trait_vector_extractor():
    extractor = TraitVectorExtractor(host="nbenthamiana")
    seq = "ATGGCTTGGTAA"
    trait = extractor.extract(seq)
    assert trait.cai_golden_set > 0.0
    assert 30.0 <= trait.global_gc_percent <= 60.0
    assert trait.homopolymer_max_run >= 1
    assert 0.0 <= trait.synthesis_penalty <= 1.0
    assert 0.0 <= trait.codon_context_prior <= 1.0


def test_discovery_slate_engine_generation():
    engine = DiscoverySlateEngine(host="nbenthamiana")
    # Small test protein: Humira LC initiation snippet (20 AA)
    target_aa = "DIQMTQSPSSLSASVGDRVT"
    slate = engine.generate_slate(
        target_aa=target_aa,
        target_name="Humira-LC-Fragment",
        mature_protein_aa_length=20,
        construct_aa_length=20,
        signal_peptide_included=False,
        novelty_class="Class_A_InDistribution",
        top_k=3,
    )

    assert slate.target_metadata.target_name == "Humira-LC-Fragment"
    assert slate.target_metadata.mature_protein_aa_length == 20
    assert slate.slate_summary.top_k_count <= 3
    assert len(slate.candidates) >= 1

    # Check candidates structure
    for cand in slate.candidates:
        assert cand.rank >= 1
        assert (
            cand.sequence_dna.startswith("G")
            or cand.sequence_dna.startswith("A")
            or cand.sequence_dna.startswith("C")
            or cand.sequence_dna.startswith("T")
        )
        assert cand.sequence_digest.startswith("sha256:")
        assert cand.trait_vector.cai_golden_set > 0.5
        assert cand.utility_score > 0.0

    # Verify distinct candidate strategies if multiple
    if len(slate.candidates) > 1:
        assert slate.slate_summary.diversity_index >= 0.0


def test_fresh_process_determinism():
    """Verify that running the slate engine twice produces byte-identical sequence digests."""
    engine1 = DiscoverySlateEngine(host="nbenthamiana")
    engine2 = DiscoverySlateEngine(host="nbenthamiana")
    target_aa = "EVQLVESGGGLVQPGRSLRL"

    slate1 = engine1.generate_slate(target_aa=target_aa, top_k=3)
    slate2 = engine2.generate_slate(target_aa=target_aa, top_k=3)

    assert len(slate1.candidates) == len(slate2.candidates)
    for c1, c2 in zip(slate1.candidates, slate2.candidates):
        assert c1.generation_contract == c2.generation_contract
        assert c1.sequence_dna == c2.sequence_dna
        assert c1.sequence_digest == c2.sequence_digest
        assert c1.utility_score == c2.utility_score
