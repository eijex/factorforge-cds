"""Unit tests for Biomolecular Structure Viewer MCP Bridge and Coordinate Mapper."""

from pathlib import Path
import pytest

from factorforge.structure.viewer_bridge import (
    build_structure_position_mapping,
    prepare_viewer_annotations,
)

SAMPLE_CANDIDATE_REGIONS = [
    {
        "start_1based": 1,
        "end_1based": 18,
        "region_type": "uniprot_signal_peptide",
        "reason": "Signal peptide cleavage",
    },
    {
        "start_1based": 47,
        "end_1based": 53,
        "region_type": "predicted_flexible_hinge",
        "reason": "EC1 loop",
    },
]

# Canonical human CD9 sequence (UniProt P21926, 228 aa)
CD9_HUMAN_228 = (
    "MPVKGGTKCIKYLLFGFNFIFWLAGIAVLAIGLWLRFDSQTKSIFEQETNNNNSSFYTGVYILIGAGALMMLVGFLGCC"
    "GAVQESQCMLGLFFGFLLVIFAIEIAAAIWGYSHKDEVIKEVQEFYKDTYNKLKTKDEPQRETLKAIHYALNCCGLAGG"
    "VEQFISDICPKKDVLETFTVKSCPDAIKEVFDNKFHIIGAVGIGIAVVMIFGMIFSMILCCAIRRNREMV"
)

# Recombinant CD9 with C-terminal 8xHis tag (236 aa, matching AlphaFold3 model)
CD9_RECOMBINANT_236 = CD9_HUMAN_228 + "HHHHHHHH"

# Arbitrary synthetic mismatch sequence
SYNTHETIC_MISMATCH_SEQ = "ACDEFGHIKLMNPQRSTVWY" * 11  # 220 aa


def test_viewer_bridge_unverified_fallback() -> None:
    """Without position mapping, annotations fall back to assumed 1-to-1 sequence indices."""
    annotations = prepare_viewer_annotations({"candidate_regions": SAMPLE_CANDIDATE_REGIONS}, auth_asym_id="A")
    assert len(annotations) == 2
    assert annotations[0]["start_auth_seq_id"] == 1
    assert annotations[0]["end_auth_seq_id"] == 18
    assert annotations[0]["coordinate_mapping"] == "ASSUMED_1_TO_1_UNVERIFIED"
    assert annotations[0]["mapping_status"] == "UNVERIFIED"
    assert "Structure Unverified" in annotations[0]["text"]


def test_viewer_bridge_with_offset_structure_mapping() -> None:
    """When structure starts at an offset (e.g. auth_seq_id 23 due to cleaved signal peptide)."""
    # 1-based sequence position 1..10 corresponds to structure auth_seq_id 23..32
    position_mapping = [
        {"protein_position": i, "auth_seq_id": i + 22, "chain_id": "A"}
        for i in range(1, 100)
    ]

    annotations = prepare_viewer_annotations(
        {"candidate_regions": SAMPLE_CANDIDATE_REGIONS},
        auth_asym_id="A",
        position_mapping=position_mapping,
    )

    # Region 1 (Seq 1..18) -> Struct 23..40
    assert annotations[0]["start_auth_seq_id"] == 23
    assert annotations[0]["end_auth_seq_id"] == 40
    assert annotations[0]["coordinate_mapping"] == "VERIFIED_STRUCTURE_AUTH_SEQ_ID"
    assert annotations[0]["mapping_status"] == "EXACT"

    # Region 2 (Seq 47..53) -> Struct 69..75
    assert annotations[1]["start_auth_seq_id"] == 69
    assert annotations[1]["end_auth_seq_id"] == 75
    assert annotations[1]["coordinate_mapping"] == "VERIFIED_STRUCTURE_AUTH_SEQ_ID"
    assert annotations[1]["mapping_status"] == "EXACT"


def test_viewer_bridge_with_missing_residues() -> None:
    """When some residues are absent in structure (e.g. disordered loop or terminal deletion)."""
    # Only residues 50..53 are present in structure
    position_mapping = [
        {"protein_position": i, "auth_seq_id": i, "chain_id": "A"}
        for i in range(50, 60)
    ]

    annotations = prepare_viewer_annotations(
        {"candidate_regions": SAMPLE_CANDIDATE_REGIONS},
        auth_asym_id="A",
        position_mapping=position_mapping,
    )

    # Region 1 (Seq 1..18) has NO residues in structure -> unmapped
    assert annotations[0]["coordinate_mapping"] == "UNMAPPED_RESIDUES_MISSING_FROM_STRUCTURE"
    assert annotations[0]["mapping_status"] == "UNMAPPED"
    assert "Missing in 3D Structure" in annotations[0]["text"]

    # Region 2 (Seq 47..53) has partial residues (50..53) -> resolved to partial present span
    assert annotations[1]["start_auth_seq_id"] == 50
    assert annotations[1]["end_auth_seq_id"] == 53
    assert annotations[1]["coordinate_mapping"] == "VERIFIED_STRUCTURE_AUTH_SEQ_ID"
    assert "PARTIAL" in annotations[1]["mapping_status"]


def test_build_structure_position_mapping_live_cif_exact_and_substring() -> None:
    """Test live MMCIF coordinate mapping against actual CD9 AlphaFold model file."""
    cif_path = Path(__file__).resolve().parents[1] / "src" / "factorforge" / "structure" / "cd9_alphafold3_model_0.cif"
    assert cif_path.exists(), f"Structure file not found: {cif_path}"

    cif_content = cif_path.read_text(encoding="utf-8")

    # 1. Native 228 aa sequence (exact substring within the 236 aa recombinant structure)
    mapping_native = build_structure_position_mapping(
        protein_sequence=CD9_HUMAN_228,
        structure_content=cif_content,
        format="cif",
        chain_id="A",
    )
    assert len(mapping_native) == 228
    assert mapping_native[0]["protein_position"] == 1
    assert mapping_native[0]["auth_seq_id"] == 1
    assert mapping_native[-1]["protein_position"] == 228
    assert mapping_native[-1]["auth_seq_id"] == 228

    # 2. Recombinant 236 aa sequence (full exact match including C-terminal His-tag)
    mapping_rec = build_structure_position_mapping(
        protein_sequence=CD9_RECOMBINANT_236,
        structure_content=cif_content,
        format="cif",
        chain_id="A",
    )
    assert len(mapping_rec) == 236
    assert mapping_rec[0]["auth_seq_id"] == 1
    assert mapping_rec[-1]["auth_seq_id"] == 236


def test_build_structure_position_mapping_rejects_synthetic_mismatch() -> None:
    """Completely different synthetic sequence must be rejected, not sequentially mapped."""
    cif_path = Path(__file__).resolve().parents[1] / "src" / "factorforge" / "structure" / "cd9_alphafold3_model_0.cif"
    assert cif_path.exists()
    cif_content = cif_path.read_text(encoding="utf-8")

    # Passing synthetic sequence must raise ValueError, NEVER return a sequential fallback
    with pytest.raises(ValueError, match="Sequence mismatch"):
        build_structure_position_mapping(
            protein_sequence=SYNTHETIC_MISMATCH_SEQ,
            structure_content=cif_content,
            format="cif",
            chain_id="A",
        )
