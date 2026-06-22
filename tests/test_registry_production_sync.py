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
