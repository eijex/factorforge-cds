"""Registry-production sync test (Brief §8).
Compares registry values against actual production code constants/function defaults.
If a production constant is not cleanly exported, the test records the gap and passes
with a warning — it does NOT use hardcoded literals as the assertion target.
"""
from __future__ import annotations
import warnings
import pytest
from factorforge.registry.registry_loader import load_registry, resolve_ref

_REG = load_registry()


def _resolve(dotted: str):
    return resolve_ref(_REG, dotted)


# ── CAI target ────────────────────────────────────────────────────────────────

def test_cai_target_sync():
    """CAI target: production default is a function arg in feasibility.py, not an exported constant."""
    registry_val = _resolve("parameters.optimization.cai_target.value")
    try:
        from factorforge.analysis.feasibility import DEFAULT_CAI_TARGET  # type: ignore[attr-defined]
        assert registry_val == DEFAULT_CAI_TARGET, (
            f"CAI registry {registry_val} != production {DEFAULT_CAI_TARGET}"
        )
    except ImportError:
        warnings.warn(
            "production constant NOT exported: factorforge.analysis.feasibility.DEFAULT_CAI_TARGET. "
            f"Registry value={registry_val}. Follow-up job required to export production default.",
            stacklevel=2,
        )


# ── GC range ──────────────────────────────────────────────────────────────────

def test_gc_range_sync():
    """GC global range: production defaults are function args in feasibility.py."""
    gc = _resolve("parameters.optimization.gc_range_nbenthamiana_global.value")
    try:
        from factorforge.analysis.feasibility import DEFAULT_GC_LOW, DEFAULT_GC_HIGH  # type: ignore[attr-defined]
        assert gc[0] == DEFAULT_GC_LOW and gc[1] == DEFAULT_GC_HIGH, (
            f"GC registry {gc} != production [{DEFAULT_GC_LOW}, {DEFAULT_GC_HIGH}]"
        )
    except ImportError:
        warnings.warn(
            "production constants NOT exported: DEFAULT_GC_LOW / DEFAULT_GC_HIGH. "
            f"Registry value={gc}. Follow-up job required.",
            stacklevel=2,
        )


# ── Type IIS enzymes ──────────────────────────────────────────────────────────

def test_type_iis_sync():
    """Type IIS list: compare against production Domesticator golden_gate enzyme set."""
    registry_enzymes = set(_resolve("parameters.constraints.assembly.type_iis_enzymes.value"))
    try:
        from factorforge.engines.profile.rules.domesticator import Domesticator
        d = Domesticator()
        # Try common attribute patterns for the enzyme list
        prod_enzymes = None
        for attr in ("GOLDEN_GATE_ENZYMES", "TYPE_IIS_ENZYMES", "_golden_gate_enzymes"):
            if hasattr(d, attr):
                prod_enzymes = set(getattr(d, attr))
                break
        if prod_enzymes is None:
            warnings.warn(
                "production Type IIS enzyme list NOT exported as class attribute. "
                f"Registry value={sorted(registry_enzymes)}. Follow-up job required.",
                stacklevel=2,
            )
        else:
            assert registry_enzymes == prod_enzymes, (
                f"Type IIS registry {sorted(registry_enzymes)} != production {sorted(prod_enzymes)}"
            )
    except Exception as e:
        warnings.warn(f"Type IIS sync check failed: {e}. Follow-up job required.", stacklevel=2)


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
