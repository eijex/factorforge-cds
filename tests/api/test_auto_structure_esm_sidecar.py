from __future__ import annotations

from unittest.mock import patch

from factorforge.structure.alphafold_server_ingest import extract_pause_positions_from_sidecar
from factorforge.structure.auto_resolver import (
    auto_resolve_structure_for_protein,
    build_4dp_sidecar_from_structure,
)
from factorforge.structure.esmfold_client import _generate_mock_pdb


def test_build_4dp_sidecar_from_pdb_extracts_pause_positions() -> None:
    protein = "MAKTNLFLFLIFSLLLSLSSAEVQLVESGG" * 2
    pdb_text = _generate_mock_pdb(protein, chain="A")
    sidecar = build_4dp_sidecar_from_structure(
        pdb_text, "pdb", chain_id="A", expected_sequence=protein
    )
    assert sidecar is not None
    assert sidecar["tunnel_offset_applied"] == 35
    # Verify compatibility with 4DP pause position extraction
    positions = extract_pause_positions_from_sidecar(sidecar, protein, allow_mock=True)
    assert len(positions) >= 1


def test_auto_resolve_structure_allow_remote_uses_esmfold_and_generates_sidecar() -> None:
    protein = "MAKTNLFLFLIFSLLLSLSSAEVQLVESGG" * 2
    mock_pdb = _generate_mock_pdb(protein, chain="A")
    with patch(
        "factorforge.structure.esmfold_client.predict_structure_esm",
        return_value=mock_pdb,
    ) as mock_predict:
        resolved = auto_resolve_structure_for_protein(protein, allow_remote=True)
    mock_predict.assert_called_once()
    assert resolved["format"] == "pdb"
    assert resolved["content"] == mock_pdb
    assert resolved["sidecar_4dp"] is not None
    assert len(resolved["ca_trace"]) == len(protein)
    assert "ESMFold" in resolved["source_label"]
