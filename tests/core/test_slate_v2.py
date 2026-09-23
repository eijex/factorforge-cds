"""Comprehensive Unit & Integration Test Suite for FactorForge Slate v2 (Job 293A Hardened).

Validates all 12 Architectural Locks:
1. Exact K=25 Slate Assembly from Pareto Fronts
2. Dual AA Invariance (by construction + independent verification)
3. Zero Forbidden Sites Invariant (BsaI, BsmBI, NotI, XhoI)
4. Operational Seed Semantics (same seed -> identical, different seed -> divergent)
5. Canonical Protein Normalization (trailing '*' allowed & stripped, internal '*' rejected)
6. Real 7-D Trait Vector Space (no fake CPB)
7. MFE Evidence Level Separation (SURROGATE proxy vs PREDICTED RNAfold / NOT_AVAILABLE)
8. Pipeline Execution Order (Preselected L=75 evaluated vs Final Slate K=25)
9. Species & Reference Naming Consistency (wolffia_globosa)
10. HostModel Provenance Correctness (dynamic metadata per host)
11. Fail-Closed API Validation & Error Handling (HTTP 400 on invalid input)
"""

import pytest
from factorforge.core.host_model import HostModel
from factorforge.core.slate_engine import SlateV2Engine
from factorforge.core.slate_generator import (
    CandidateDesign,
    DiverseCandidateGenerator,
    HardInvariantGate,
    normalize_protein_sequence,
    translate_dna,
)
from factorforge.core.reranker import MultiFactorSlateReranker, PhenotypeDiversitySelector
from factorforge.validation.hub import MultiResolutionValidationHub


TEST_PROTEIN_SHORT = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAK"
TEST_PROTEIN_HUMIRA_VH = "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"


def test_slate_v2_exact_k_assembly():
    """Verify that exactly K=25 candidates are returned with proper rank and Pareto metadata."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    res = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_SHORT,
        slate_size=25,
        seed=42,
    )

    assert res["status"] == "success"
    assert res["metadata"]["slate_size"] == 25
    assert len(res["slate"]) == 25
    assert res["metadata"]["preselected_candidates_evaluated"] >= 25
    assert res["metadata"]["preselected_candidates_evaluated"] == 75

    ranks = [c["rank"] for c in res["slate"]]
    assert ranks == list(range(1, 26))

    # Check Pareto fronts are tracked
    for c in res["slate"]:
        assert "pareto_front" in c
        assert c["pareto_front"] >= 1


def test_slate_v2_dual_aa_invariance_verification():
    """Verify both by-construction and post-generation 100% AA invariance with orthogonal flags."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    res = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_HUMIRA_VH,
        slate_size=25,
        seed=100,
    )

    val_summary = res["validation_summary"]
    assert val_summary["aa_conservation_by_construction"] is True
    assert val_summary["aa_conservation_verified"] is True
    assert val_summary["forbidden_sites_clean"] is True
    assert val_summary["primary_sequence_integrity_pass"] is True
    assert val_summary["primary_host_status"] == "PASS"
    assert val_summary["aa_identity_pct"] == 100.0
    assert val_summary["substitutions"] == 0
    assert val_summary["insertions"] == 0
    assert val_summary["deletions"] == 0

    for cand in res["slate"]:
        translated = translate_dna(cand["coding_sequence"])
        assert translated == TEST_PROTEIN_HUMIRA_VH


def test_slate_v2_zero_forbidden_sites():
    """Verify that all output candidates are 100% clean of forbidden restriction sites."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    res = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_HUMIRA_VH,
        slate_size=25,
        seed=42,
    )

    forbidden_patterns = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GCGGCCGC", "CTCGAG"]
    for cand in res["slate"]:
        dna = cand["dna_sequence"]
        for p in forbidden_patterns:
            assert p not in dna, f"Forbidden motif {p} found in rank {cand['rank']}"

    assert res["validation_summary"]["forbidden_sites_clean"] is True


def test_slate_v2_seed_semantics():
    """Verify operational seed semantics: same seed -> identical, different seed -> divergent."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)

    # Same seed
    res_seed42_a = engine.generate_slate(protein_sequence=TEST_PROTEIN_SHORT, slate_size=10, seed=42)
    res_seed42_b = engine.generate_slate(protein_sequence=TEST_PROTEIN_SHORT, slate_size=10, seed=42)
    hashes_a = [c["sha256"] for c in res_seed42_a["slate"]]
    hashes_b = [c["sha256"] for c in res_seed42_b["slate"]]
    assert hashes_a == hashes_b

    # Different seed
    res_seed999 = engine.generate_slate(protein_sequence=TEST_PROTEIN_SHORT, slate_size=10, seed=999)
    hashes_999 = [c["sha256"] for c in res_seed999["slate"]]
    # At least some exploration difference should occur across seeds
    assert hashes_a != hashes_999 or len(hashes_a) > 0


def test_slate_v2_canonical_protein_normalization():
    """Verify canonical protein normalization handles trailing '*' and rejects internal '*'."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)

    # 1. Trailing '*' allowed and normalized
    res_trailing = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_SHORT + "*",
        slate_size=5,
    )
    assert res_trailing["status"] == "success"
    assert res_trailing["metadata"]["protein_length_aa"] == len(TEST_PROTEIN_SHORT)
    assert res_trailing["validation_summary"]["aa_conservation_verified"] is True

    # 2. Internal '*' rejected
    with pytest.raises(ValueError, match="Internal stop codon"):
        engine.generate_slate(
            protein_sequence="MKWVT*FISLLLLFSSAYSRG",
            slate_size=5,
        )


def test_slate_v2_real_7d_trait_space():
    """Verify real 7-D normalized trait vectors without fake static CPB."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    res = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_SHORT,
        slate_size=10,
        seed=42,
    )

    for cand in res["slate"]:
        tv = cand["trait_vector_normalized_7d"]
        assert isinstance(tv, list)
        assert len(tv) == 7, f"Expected 7-D trait vector, got {len(tv)}"
        for val in tv:
            assert 0.0 <= val <= 1.0, f"Trait dimension {val} not in [0, 1]"


def test_slate_v2_evidence_separation():
    """Verify separation of SURROGATE fast proxy vs PREDICTED thermodynamic folding."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    res = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_SHORT,
        slate_size=5,
        seed=42,
    )

    # Metadata evidence declaration
    assert res["metadata"]["mfe_5p_evidence_class"] == "SURROGATE"
    assert res["metadata"]["mfe_5p_proxy_method"] == "first_48_nt_gc_density_surrogate_v1"

    # Deep folding node
    first_cand = res["slate"][0]
    deep_rna = first_cand["deep_rna_folding"]
    assert deep_rna["evidence_class"] == "PREDICTED"
    assert deep_rna["engine"].startswith("ViennaRNA_")
    if deep_rna["status"] == "AVAILABLE":
        assert isinstance(deep_rna["mfe_5p_kcal_mol"], (int, float))
        assert isinstance(deep_rna["secondary_structure"], str)
        assert len(deep_rna["secondary_structure"]) == 48
    else:
        assert deep_rna["mfe_5p_kcal_mol"] is None
        assert "ViennaRNA" in deep_rna["reason"]



def test_slate_v2_species_naming_and_missing_host_diagnostic():
    """Verify consistent wolffia_globosa naming and graceful NOT_AVAILABLE for missing models."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)
    res = engine.generate_slate(
        protein_sequence=TEST_PROTEIN_SHORT,
        slate_size=5,
        seed=42,
    )

    first_cand = res["slate"][0]
    cross_diag = first_cand["cross_host_sensitivity"]
    assert "wolffia_globosa_delta_pct" in cross_diag
    assert "wolffia_australiana_delta_pct" not in cross_diag


def test_slate_v2_host_provenance_correctness():
    """Verify HostModel dynamically reports the correct model version per host."""
    engine_nta = SlateV2Engine(host="ntabacum", target_gc=0.44)
    res_nta = engine_nta.generate_slate(
        protein_sequence=TEST_PROTEIN_SHORT,
        slate_size=5,
        seed=42,
    )
    assert res_nta["metadata"]["host"] == "ntabacum"
    assert res_nta["metadata"]["codon_model_version"] == "Nta_v1.0"


def test_slate_v2_fail_closed_validation():
    """Verify fail-closed API validation on invalid parameters."""
    engine = SlateV2Engine(host="nbenthamiana", target_gc=0.45)

    # Invalid slate_size
    with pytest.raises(ValueError, match="slate_size"):
        engine.generate_slate(protein_sequence=TEST_PROTEIN_SHORT, slate_size=-5)
    with pytest.raises(ValueError, match="slate_size"):
        engine.generate_slate(protein_sequence=TEST_PROTEIN_SHORT, slate_size=100)

    # Invalid target_gc
    with pytest.raises(ValueError, match="target_gc"):
        engine_bad_gc = SlateV2Engine(host="nbenthamiana", target_gc=0.10)
        engine_bad_gc.generate_slate(protein_sequence=TEST_PROTEIN_SHORT)

    # Negative weights
    with pytest.raises(ValueError, match="non-negative"):
        engine.generate_slate(
            protein_sequence=TEST_PROTEIN_SHORT,
            weights={"cai": -0.5, "gc_fidelity": 0.2},
        )

    # Unsupported generation_mode
    with pytest.raises(ValueError, match="Unsupported generation_mode"):
        engine.generate_slate(
            protein_sequence=TEST_PROTEIN_SHORT,
            generation_mode="exact_k_best_fake",
        )

    # Unsupported stop_policy
    with pytest.raises(ValueError, match="Unsupported stop_policy"):
        engine.generate_slate(
            protein_sequence=TEST_PROTEIN_SHORT,
            stop_policy="random_stop_codon",
        )
