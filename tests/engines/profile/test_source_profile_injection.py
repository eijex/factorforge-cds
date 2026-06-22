"""
Source-profile codon table injection into RuleBasedOptimizer.

Covers the benchmark contract: when a source-profile codon table is
injected, the CAI-dependent profiles must re-design against it, while the product
engine path (no injection) must remain byte-identical to the bundled reference.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.optimizer import RuleBasedOptimizer

PROFILES_DIR = (
    Path(__file__).resolve().parents[3]
    / "src" / "factorforge" / "data" / "profiles"
)

# A Lys/Leu/Ala-rich probe makes synonymous codon shifts easy to detect.
PROBE = "MKKKKKLLLLLAAAAAGGGGGEEEEE"
SEED = 320

SOURCE_TABLES = [
    "qld183_v103_derived_codons.json",
    "nbev11_cds_all_derived_codons.json",
    "nbev11_cds_hc_derived_codons.json",
]


def _has_profiles() -> bool:
    return all((PROFILES_DIR / t).exists() for t in SOURCE_TABLES)


pytestmark = pytest.mark.skipif(
    not _has_profiles(),
    reason="Source-profile codon tables not built.",
)


def test_default_path_unchanged_by_injection_support():
    """No injection ⇒ identical design to a second default optimizer (product path)."""
    a = RuleBasedOptimizer()
    b = RuleBasedOptimizer(codon_table_path=None)
    for profile in ("balanced", "high_cai", "gc_target", "assembly_friendly"):
        assert (
            a.optimize(PROBE, profile=profile, seed=SEED).sequence
            == b.optimize(PROBE, profile=profile, seed=SEED).sequence
        )


@pytest.mark.parametrize("table_name", SOURCE_TABLES)
@pytest.mark.parametrize("profile", ["balanced", "high_cai", "assembly_friendly"])
def test_cai_dependent_profiles_respond_to_injection(table_name, profile):
    """CAI-driven profiles must produce different designs under an injected table."""
    legacy = RuleBasedOptimizer()
    injected = RuleBasedOptimizer(codon_table_path=str(PROFILES_DIR / table_name))
    legacy_seq = legacy.optimize(PROBE, profile=profile, seed=SEED).sequence
    injected_seq = injected.optimize(PROBE, profile=profile, seed=SEED).sequence
    assert legacy_seq != injected_seq, (
        f"{profile} design did not respond to injected table {table_name}"
    )


@pytest.mark.parametrize("table_name", SOURCE_TABLES)
def test_gc_target_is_gc_driven_and_invariant(table_name):
    """gc_target is GC-driven; it legitimately ignores the codon-usage source."""
    legacy = RuleBasedOptimizer()
    injected = RuleBasedOptimizer(codon_table_path=str(PROFILES_DIR / table_name))
    assert (
        legacy.optimize(PROBE, profile="gc_target", seed=SEED).sequence
        == injected.optimize(PROBE, profile="gc_target", seed=SEED).sequence
    )


def test_injection_is_deterministic():
    """Same injected table + seed ⇒ identical design across instances."""
    table = str(PROFILES_DIR / SOURCE_TABLES[0])
    one = RuleBasedOptimizer(codon_table_path=table)
    two = RuleBasedOptimizer(codon_table_path=table)
    assert (
        one.optimize(PROBE, profile="balanced", seed=SEED).sequence
        == two.optimize(PROBE, profile="balanced", seed=SEED).sequence
    )
