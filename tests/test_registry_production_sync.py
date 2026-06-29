"""Registry-production sync test (Brief §8).
Compares registry values against actual production code constants/function defaults.
Production constants are exported as named module/class attributes,
so these checks are strict assertions against the registry (single source of truth).
"""
from __future__ import annotations
from factorforge.registry.registry_loader import load_registry, resolve_ref

_REG = load_registry()


def _resolve(dotted: str):
    return resolve_ref(_REG, dotted)


# ── CAI target ────────────────────────────────────────────────────────────────

def test_cai_target_sync():
    """CAI target: registry value must match the exported production default."""
    registry_val = _resolve("parameters.optimization.cai_target.value")
    from factorforge.analysis.feasibility import DEFAULT_CAI_TARGET
    assert registry_val == DEFAULT_CAI_TARGET, (
        f"CAI registry {registry_val} != production {DEFAULT_CAI_TARGET}"
    )


# ── GC range ──────────────────────────────────────────────────────────────────

def test_gc_range_sync():
    """GC global range: registry value must match the exported production defaults."""
    gc = _resolve("parameters.optimization.gc_range_nbenthamiana_global.value")
    from factorforge.analysis.feasibility import DEFAULT_GC_LOW, DEFAULT_GC_HIGH
    assert gc[0] == DEFAULT_GC_LOW and gc[1] == DEFAULT_GC_HIGH, (
        f"GC registry {gc} != production [{DEFAULT_GC_LOW}, {DEFAULT_GC_HIGH}]"
    )


# ── Type IIS enzymes ──────────────────────────────────────────────────────────

def test_type_iis_sync():
    """Type IIS list: registry value must match the exported Domesticator enzyme set."""
    registry_enzymes = set(_resolve("parameters.constraints.assembly.type_iis_enzymes.value"))
    from factorforge.engines.profile.rules.domesticator import Domesticator
    prod_enzymes = set(Domesticator.GOLDEN_GATE_ENZYMES)
    assert registry_enzymes == prod_enzymes, (
        f"Type IIS registry {sorted(registry_enzymes)} != production {sorted(prod_enzymes)}"
    )


# ── BbsI exclusion ───────────────────────────────────────────────────────────

def test_bbsi_not_in_canonical_enzyme_set():
    """BbsI must not appear in the active Type IIS enzyme list.

    BpiI is the canonical FactorForge label for the GAAGAC target.
    BbsI is an isoschizomer alias and must never re-enter as an active scan
    target — it was explicitly removed to prevent duplicate scanning.
    """
    from factorforge.engines.profile.rules.domesticator import Domesticator
    assert "BbsI" not in set(Domesticator.GOLDEN_GATE_ENZYMES), (
        "BbsI must not be in GOLDEN_GATE_ENZYMES — use BpiI (canonical label)"
    )


# ── AA identity / internal stop ───────────────────────────────────────────────

def test_aa_identity_policy_sync():
    """AA identity: validate_cds_output must require exact match (1.0)."""
    registry_val = _resolve("parameters.constraints.biological.aa_identity_target.value")
    from factorforge.utils.sequence_validator import validate_cds_output
    # A mismatched sequence should fail
    result = validate_cds_output("MK", "ATGAAAACC")  # ATG=M AAA=K ACC=T → MKT ≠ MK
    assert not result["passed"], "validator should fail on AA mismatch"
    # registry says identity must be 1.0 — passing case confirms
    ok = validate_cds_output("MKT", "ATGAAAACC")
    assert ok["passed"] and ok["aa_identity"] == registry_val


# ── codon_reference source-of-truth sync (Job 168 / v3.3.0, _analysis/025) ────

def test_codon_reference_active_sync_with_active_reference_file():
    """registry's codon_reference.active block must match
    data/reference/active_codon_reference.json (the file run_benchmark.py
    reads at runtime) — both describe "what is the current production
    default", and must not drift apart."""
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    registry_active = _resolve("codon_reference.active")
    active_ref = json.loads(
        (root / "data" / "reference" / "active_codon_reference.json").read_text(encoding="utf-8")
    )
    assert registry_active["id"] == active_ref["active_codon_table_id"]
    assert registry_active["asset_type"] == active_ref["active_asset_type"]
    assert registry_active["codon_reference_contract_version"] == active_ref[
        "codon_reference_contract_version"
    ]
    # active_codon_reference.json doesn't carry sha256 directly — cross-check
    # via the file it points at instead.
    import hashlib
    table_path = root / active_ref["file"]
    assert registry_active["sha256"] == hashlib.sha256(table_path.read_bytes()).hexdigest()


def test_codon_reference_active_sync_with_legacy_manifest():
    """registry's codon_reference.active block must match the schema-conformant
    legacy (v1) manifest file's facts (asset_type, sha256, source_status) — v1
    is the production default again as of the v3.2.7 GC-band/codon-reference
    revert, pending an MFE re-sensitivity + 2x2 factorial recheck of v2."""
    import json
    from pathlib import Path

    registry_active = _resolve("codon_reference.active")
    manifest = json.loads(
        (Path(__file__).resolve().parents[1] / "data" / "reference"
         / "codon_table_manifest.json").read_text(encoding="utf-8")
    )
    assert registry_active["id"] == manifest["codon_table_id"]
    assert registry_active["sha256"] == manifest["sha256"]
    assert registry_active["asset_type"] == manifest["asset_type"]
    assert registry_active["source_status"] == manifest["source_status"]


def test_codon_reference_candidate_sync_with_v2_manifest():
    """registry's codon_reference.candidate block (the provisionally
    un-promoted v2 asset) must stay in sync with its own manifest file even
    while not active."""
    import json
    from pathlib import Path

    registry_candidate = _resolve("codon_reference.candidate")
    manifest = json.loads(
        (Path(__file__).resolve().parents[1] / "data" / "reference"
         / "codon_table_manifest_nbev11_hc_v2.json").read_text(encoding="utf-8")
    )
    assert registry_candidate["id"] == manifest["codon_table_id"]
    assert registry_candidate["sha256"] == manifest["sha256"]
    assert registry_candidate["asset_type"] == manifest["asset_type"]
    assert registry_candidate["source_status"] == manifest["source_status"]


def test_codon_reference_active_table_sha256_matches_production_default():
    """The sha256 recorded for the active codon_reference must match the
    actual file the production engine resolves to by default."""
    from factorforge.engines.profile.utils import get_data_path, resolve_host_codon_table_path
    import hashlib

    registry_active = _resolve("codon_reference.active")
    resolved_path = resolve_host_codon_table_path("nbenthamiana", get_data_path())
    actual_sha256 = hashlib.sha256(resolved_path.read_bytes()).hexdigest()
    assert registry_active["sha256"] == actual_sha256
