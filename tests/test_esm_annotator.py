"""Tests for ESM-2 Neural Contact & Flexibility Propensity Adapter."""

import pytest
from factorforge.structure.esm_annotator import (
    ESM_ANNOTATOR_VERSION,
    ESM_LIMITATIONS_STATEMENT,
    ESMStructureAnnotator,
    recommend_esm_pause_regions,
)

ADALIMUMAB_VH = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
)


def test_esm_annotator_live_inference() -> None:
    annotator = ESMStructureAnnotator(model_name="esm2_t6_8M_UR50D")
    if not annotator.is_available():
        pytest.skip("fair-esm or weights not available in environment")

    report = annotator.annotate(ADALIMUMAB_VH)

    # Invariants
    assert report.status == "AVAILABLE"
    assert report.protein_length == len(ADALIMUMAB_VH)
    assert report.reviewed is False  # Mandatory human-in-the-loop invariant
    assert report.tool_version == ESM_ANNOTATOR_VERSION
    assert report.limitations == ESM_LIMITATIONS_STATEMENT
    assert report.provenance_basis == "UNVERIFIED_PUBLIC_CORPUS"
    assert report.biological_validation == "NOT_AVAILABLE"

    # Contact metrics
    assert report.contact_density_min > 0.0
    assert report.contact_density_max > report.contact_density_min
    assert report.contact_density_mean > 0.0

    # Regions & positions
    assert len(report.candidate_regions) >= 1
    assert len(report.recommended_positions) >= 1
    assert all(1 <= p <= len(ADALIMUMAB_VH) for p in report.recommended_positions)

    for region in report.candidate_regions:
        assert 1 <= region.start_1based <= region.end_1based <= len(ADALIMUMAB_VH)
        assert region.length == region.end_1based - region.start_1based + 1
        assert region.region_type == "predicted_flexible_hinge"


def test_esm_annotator_empty_sequence() -> None:
    annotator = ESMStructureAnnotator()
    report = annotator.annotate("")
    assert report.status == "ERROR"
    assert report.protein_length == 0
    assert "empty" in (report.reason or "").lower()


def test_esm_annotator_convenience_function() -> None:
    annotator = ESMStructureAnnotator()
    if not annotator.is_available():
        pytest.skip("fair-esm not available")

    # Short test peptide
    pep = "ACDEFGHIKLMNPQRSTVWYACDEF"
    report = recommend_esm_pause_regions(pep)
    assert report.status == "AVAILABLE"
    assert report.protein_length == 25
    assert report.reviewed is False


def test_esm_annotator_missing_model_graceful_fallback() -> None:
    annotator = ESMStructureAnnotator(model_name="non_existent_esm_model_xyz")
    assert annotator.is_available() is False
    report = annotator.annotate("ACDEF")
    assert report.status == "UNAVAILABLE"
    assert "could not be loaded" in (report.reason or "")
