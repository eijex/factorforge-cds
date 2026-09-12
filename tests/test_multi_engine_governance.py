"""Multi-Engine Governance & Provenance Test Suite for FactorForge v3.5.0."""

import pytest
from factorforge.engines.dp_v2 import DPV2Optimizer
from factorforge.registry.versioning import (
    product_version,
    engine_version,
    engine_generation,
    engine_status,
    engine_runtime_state,
    public_version_metadata,
)
from factorforge.analysis.metrics import load_codon_usage_table


def test_version_manifest_governance():
    """Verify product and multi-generation engine SemVer decoupling."""
    assert product_version() == "3.5.0"
    
    # Gen 1 (Rule)
    assert engine_generation("profile") == 1
    assert engine_version("profile") == "1.0.0"
    assert engine_status("profile") == "stable"
    
    # Gen 2 (DP v2)
    assert engine_generation("dp") == 2
    assert engine_version("dp") == "2.0.1"
    assert engine_status("dp") == "stable"
    
    # Gen 3 (sLLM Hybrid)
    assert engine_generation("slm") == 3
    assert engine_version("slm") == "0.1.0-preview.1"
    assert engine_status("slm") == "research_preview"
    assert engine_runtime_state("slm") == "scaffold"


def test_dp_v2_deterministic_tie_breaking():
    """Verify 100% byte-identical repeatable output across multiple runs."""
    peptide = "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSG"
    table = load_codon_usage_table()
    optimizer = DPV2Optimizer()
    
    res1 = optimizer.optimize(peptide, codon_weights=table.codon_weights, target_gc_min=0.40, target_gc_max=0.47)
    res2 = optimizer.optimize(peptide, codon_weights=table.codon_weights, target_gc_min=0.40, target_gc_max=0.47)
    res3 = optimizer.optimize(peptide, codon_weights=table.codon_weights, target_gc_min=0.40, target_gc_max=0.47)
    
    assert res1["sequence"] == res2["sequence"] == res3["sequence"]
    assert res1["cai"] == res2["cai"] == res3["cai"]
    assert res1["gc_percent"] == res2["gc_percent"] == res3["gc_percent"]


def test_api_provenance_envelope():
    """Verify /api/optimize response structure embeds complete engine provenance."""
    from api.optimize import handler
    
    instance = handler.__new__(handler)
    peptide = "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIY"
    
    response = instance.optimize_sequence(
        sequence=peptide,
        profile="balanced",
        use_template=False,
        kozak=False,
        dinuc=False,
        objective="feasibility_best",
        host_profile="Nicotiana benthamiana",
        host="nbenthamiana",
        return_candidates=True,
        constraints={"gc_min": 40.0, "gc_max": 47.0, "cai_target": 0.85},
    )
    
    assert response["success"] is True
    assert "provenance" in response
    prov = response["provenance"]
    assert prov["product_version"] == "3.5.0"
    assert prov["engine_id"] == "dp"
    assert prov["engine_generation"] == 2
    assert prov["engine_version"] == "2.0.1"
    assert prov["engine_status"] == "stable"
    assert "run_id" in prov
    assert "input_sequence_hash" in prov
    assert "output_cds_hash" in prov