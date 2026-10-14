"""Tests for 3D Structure Pause/Linker Recommendation Adapter (Job 395)."""

import pytest
from Bio.SeqUtils import seq3

from factorforge.structure.auto_annotator import (
    LIMITATIONS_STATEMENT,
    RECOMMENDER_VERSION,
    recommend_structure_pause_regions,
)


def _generate_synthetic_pdb(
    sequence: str,
    b_factors: list[float] | None = None,
    chain: str = "A",
    header: str = "HEADER    SYNTHETIC STRUCTURE",
) -> str:
    """Generate a minimal valid PDB with C-alpha coordinates."""
    lines = [header]
    for i, aa in enumerate(sequence, start=1):
        res3 = seq3(aa).upper()
        bf = b_factors[i - 1] if b_factors else 30.0
        line = f"ATOM  {i:5d}  CA  {res3:3s} {chain}{i:4d}    {i*3.8:8.3f}   0.000   0.000  1.00 {bf:6.2f}           C"
        lines.append(line)
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines)


def test_auto_annotator_exact_match_pass() -> None:
    # 25 amino acid synthetic sequence
    protein = "ACDEFGHIKLMNPQRSTVWYACDEF"
    # Residues 10-15 have elevated B-factors (flexible loop)
    b_factors = [20.0] * 9 + [65.0] * 6 + [20.0] * 10
    pdb_text = _generate_synthetic_pdb(protein, b_factors)

    report = recommend_structure_pause_regions(
        protein=protein,
        structure_content=pdb_text,
        format="pdb",
        chain_id="A",
        is_alphafold=False,
    )

    assert report.exact_sequence_match is True
    assert report.protein_length == 25
    assert report.reviewed is False  # Mandatory human review invariant
    assert report.tool_version == RECOMMENDER_VERSION
    assert report.limitations == LIMITATIONS_STATEMENT
    assert len(report.candidate_regions) >= 1
    assert all(1 <= p <= 25 for p in report.recommended_positions)
    assert report.coverage > 0.0


def test_auto_annotator_sequence_mismatch_raises_error() -> None:
    protein = "ACDEFGHIKL"
    wrong_protein = "ACDEFGHIKM"  # L -> M mutation
    pdb_text = _generate_synthetic_pdb(wrong_protein)

    with pytest.raises(ValueError) as exc:
        recommend_structure_pause_regions(
            protein=protein,
            structure_content=pdb_text,
            format="pdb",
            chain_id="A",
        )
    assert "does not match input protein" in str(exc.value)
    assert "Exact full-chain match is required" in str(exc.value)


def test_auto_annotator_alphafold_provenance_detection() -> None:
    protein = "ACDEFGHIKLMNPQRSTVWY"
    # Low pLDDT (<70) at positions 5..10
    plddt = [90.0] * 4 + [45.0] * 6 + [90.0] * 10
    af_pdb_text = _generate_synthetic_pdb(
        protein, plddt, header="HEADER    AlphaFold model prediction AF-TEST-F1"
    )

    report = recommend_structure_pause_regions(
        protein=protein,
        structure_content=af_pdb_text,
        format="pdb",
    )

    assert report.is_alphafold_provenance is True
    assert "AlphaFold per-residue pLDDT" in report.confidence_basis
    assert report.reviewed is False

    # Check candidate region detected
    low_plddt_regions = [r for r in report.candidate_regions if r.start_1based == 5 and r.end_1based == 10]
    assert len(low_plddt_regions) == 1
    assert low_plddt_regions[0].length == 6


def test_auto_annotator_experimental_pdb_provenance() -> None:
    protein = "ACDEFGHIKLMNPQRSTVWY"
    exp_pdb_text = _generate_synthetic_pdb(
        protein, header="HEADER    CRYSTAL STRUCTURE 1XYZ"
    )

    report = recommend_structure_pause_regions(
        protein=protein,
        structure_content=exp_pdb_text,
        format="pdb",
        is_alphafold=False,
    )

    assert report.is_alphafold_provenance is False
    assert "Experimental crystallographic/NMR B-factors" in report.confidence_basis
    assert report.reviewed is False
