"""Biological regression tests for the FactorForge v3 fallback path."""

import pytest

from factorforge.engines.v3.metrics import compute_cai, compute_gc, load_codon_usage_table
from factorforge.engines.v3.pipeline import V3Pipeline

GFP = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSR"
    "YPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYN"
    "SHVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVL"
    "LEFVTAAGITHGMDELYK"
)

CD47 = (
    "MWPLVAALLLGSACCGSAQLLFNKTKSVEHSDGDLVNEVDGSNFTVSLEPGGRRITMQLKPKDGEFIQSPTR"
    "TLDQFTFVQLNESKEVEGMAYRMV"
)

POLYA_PATTERNS = {"AATAAA", "ATTAAA", "AGTAAA", "AATACA", "AAGAAA", "AATGAA"}


@pytest.fixture(scope="module")
def pipeline() -> V3Pipeline:
    return V3Pipeline()


@pytest.fixture(scope="module")
def codon_usage():
    return load_codon_usage_table()


class TestGFP:
    def test_output_length_matches_protein(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(GFP)
        assert len(result.sequence) == len(GFP) * 3

    def test_cai_minimum(self, pipeline: V3Pipeline, codon_usage) -> None:
        result = pipeline.run(GFP)
        cai = compute_cai(result.sequence, codon_usage)
        assert cai >= 0.80, f"GFP CAI {cai:.3f} < 0.80"

    def test_gc_content_range(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(GFP)
        gc = compute_gc(result.sequence)
        assert 41.0 <= gc <= 44.0, f"GFP GC {gc:.1f}% outside 41-44% range"

    def test_no_polya_signals(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(GFP)
        found = [pattern for pattern in POLYA_PATTERNS if pattern in result.sequence]
        assert not found, f"PolyA signals found in GFP output: {found}"

    def test_valid_dna_sequence(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(GFP)
        assert all(base in "ATGC" for base in result.sequence)
        assert len(result.sequence) % 3 == 0


class TestCD47:
    def test_output_length_matches_protein(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(CD47)
        assert len(result.sequence) == len(CD47) * 3

    def test_cai_minimum(self, pipeline: V3Pipeline, codon_usage) -> None:
        result = pipeline.run(CD47)
        cai = compute_cai(result.sequence, codon_usage)
        assert cai >= 0.70, f"CD47 CAI {cai:.3f} < 0.70"

    def test_gc_content_range(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(CD47)
        gc = compute_gc(result.sequence)
        assert 41.0 <= gc <= 44.0, f"CD47 GC {gc:.1f}% outside 41-44% range"

    def test_no_polya_signals(self, pipeline: V3Pipeline) -> None:
        result = pipeline.run(CD47)
        found = [pattern for pattern in POLYA_PATTERNS if pattern in result.sequence]
        assert not found, f"PolyA signals found in CD47 output: {found}"
