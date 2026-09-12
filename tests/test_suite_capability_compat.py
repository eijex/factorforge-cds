# factorforge/tests/test_suite_capability_compat.py
"""Tests for Suite Capability Compatibility & 3-Tier Manifest Integrity."""

import json
from pathlib import Path
import pytest
from factorforge.profiles.profile_snapshot import HostProfileSnapshot


def test_host_profile_snapshot_integrity() -> None:
    sample_data = {
        "organism": "Nicotiana benthamiana",
        "snapshot_version": "nbent-2026.09",
        "description": "PlantForm Certified N. benthamiana Codon Bias Table",
        "codon_table": {"TTT": 0.85, "TTC": 0.15, "ATG": 1.0},
        "frequency_table": {"TTT": 0.024, "TTC": 0.018, "ATG": 0.021},
        "metadata": {"source": "PlantForm-RNA-seq-2026-08"},
    }
    
    # Calculate canonical digest
    digest = HostProfileSnapshot.compute_profile_digest(
        organism=sample_data["organism"],
        snapshot_version=sample_data["snapshot_version"],
        codon_table=sample_data["codon_table"],
        frequency_table=sample_data["frequency_table"],
    )
    sample_data["digest"] = digest
    
    # Load snapshot
    snapshot = HostProfileSnapshot.from_dict(sample_data)
    assert snapshot.digest == digest
    assert snapshot.organism == "Nicotiana benthamiana"
    
    # Tampering with codon table should fail fast
    tampered_data = dict(sample_data)
    tampered_data["codon_table"] = {"TTT": 0.99, "TTC": 0.01, "ATG": 1.0}
    with pytest.raises(ValueError, match="integrity mismatch"):
        HostProfileSnapshot.from_dict(tampered_data)


def test_suite_release_lock_schema_structure() -> None:
    schema_path = Path("c:/Work/eijex/eijex-workspace/manifests/schemas/suite_release_lock.schema.json")
    assert schema_path.exists()
    
    with open(schema_path, "r", encoding="utf-8") as f:
        schema = json.load(f)
        
    assert "components" in schema["required"]
    assert "scoring_specification" in schema["required"]
    assert "requires" in schema["required"]


def test_run_manifest_schema_structure() -> None:
    schema_path = Path("c:/Work/eijex/eijex-workspace/manifests/schemas/run_manifest.schema.json")
    assert schema_path.exists()
    
    with open(schema_path, "r", encoding="utf-8") as f:
        schema = json.load(f)
        
    assert "resolved_components" in schema["required"]
    assert "environment" in schema["required"]
    assert "output" in schema["required"]
