"""Tests for 5' Translation Initiation Window MFE Observation (Job 395)."""

import pytest

from factorforge.guards.initiation_mfe import (
    DEFAULT_INITIATION_WINDOW_NT,
    DEFAULT_TEMPERATURE_C,
    evaluate_initiation_mfe,
)


def test_initiation_mfe_scope_and_coordinates_without_flank() -> None:
    cds = "ATG" + "GCA" * 20  # 63 nt
    res = evaluate_initiation_mfe(cds)
    
    assert res.input_scope == "cds_5p_only"
    assert res.window_start == 0
    assert res.window_end == DEFAULT_INITIATION_WINDOW_NT
    assert res.folded_sequence_length == DEFAULT_INITIATION_WINDOW_NT
    assert res.temperature_c == DEFAULT_TEMPERATURE_C
    assert "Thermodynamic MFE is an in-silico RNA secondary structure proxy" in res.truthfulness_note
    
    # Truthful status: either AVAILABLE with float or NOT_AVAILABLE with None
    if res.status == "AVAILABLE":
        assert isinstance(res.mfe_kcal_mol, float)
        assert res.secondary_structure is not None
    else:
        assert res.status == "NOT_AVAILABLE"
        assert res.mfe_kcal_mol is None
        assert res.secondary_structure is None
        assert "not installed" in res.reason


def test_initiation_mfe_scope_and_coordinates_with_flank() -> None:
    flank = "GGTC"  # 4 nt
    cds = "ATG" + "GCA" * 20
    res = evaluate_initiation_mfe(cds, upstream_context=flank)
    
    assert res.input_scope == "upstream_flank_plus_cds_5p"
    assert res.window_start == -4
    assert res.window_end == DEFAULT_INITIATION_WINDOW_NT - 4
    assert res.folded_sequence_length == DEFAULT_INITIATION_WINDOW_NT
    
    # Must never invent a fake -4 Kozak if not provided, but since provided, coordinates match


def test_initiation_mfe_short_cds_handling() -> None:
    short_cds = "ATGGCATGA"  # 9 nt
    res = evaluate_initiation_mfe(short_cds)
    
    assert res.input_scope == "cds_5p_only"
    assert res.window_start == 0
    assert res.window_end == 9
    assert res.folded_sequence_length == 9
