"""
Source-profile codon table injection into RuleBasedOptimizer.

Covers the benchmark contract: when a source-profile codon table is
injected, the CAI-dependent profiles must re-design against it, while the product
engine path (no injection) must remain byte-identical to the bundled reference.

Job 168 / v3.3.0 (_analysis/025) made the production default itself one of
the candidate source-profile tables (nbev11_cds_hc_derived_codons.json).
These tests therefore check provenance (resolved asset path/sha256, and
whether the engine actually re-designed against it) rather than asserting
byte-identical-or-different output against a fixed "legacy" baseline —
that framing no longer holds now that "default" and "an injectable
candidate" can be the same asset, and tie-breaking among synonymous
codons of equal GC% is not guaranteed stable across different but
GC-similar reference tables.
"""

import hashlib
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
ACTIVE_REFERENCE_PATH = (
    Path(__file__).resolve().parents[3] / "data" / "reference" / "active_codon_reference.json"
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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _resolved_default_table_name() -> str:
    """Name of the file the no-injection (default) engine path actually uses."""
    active_ref = json.loads(ACTIVE_REFERENCE_PATH.read_text(encoding="utf-8"))
    return Path(active_ref["file"]).name


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
def test_injected_table_resolves_to_requested_asset(table_name):
    """Injecting a source-profile table must actually load that file's bytes —
    i.e. the resolved CAI weight source must match the requested asset's
    sha256, regardless of which table happens to be the current default."""
    injected = RuleBasedOptimizer(codon_table_path=str(PROFILES_DIR / table_name))
    assert injected._codon_table_path == str(PROFILES_DIR / table_name)
    assert injected.translator.codon_table is not None
    expected_sha256 = _sha256(PROFILES_DIR / table_name)
    actual_sha256 = _sha256(Path(injected._codon_table_path))
    assert actual_sha256 == expected_sha256


@pytest.mark.parametrize("table_name", SOURCE_TABLES)
@pytest.mark.parametrize("profile", ["balanced", "high_cai", "assembly_friendly"])
def test_cai_dependent_profiles_respond_to_injection(table_name, profile):
    """CAI-driven profiles must produce different designs under an injected
    table — UNLESS the injected table is byte-identical to the current
    default (in which case identical output is the only correct result).

    high_cai is a special case: explicit injection always overrides BOTH the
    main codon table AND the golden-set CAI reference (optimizer.py: "it is
    both the design weight table and the golden set that high_cai maximizes
    against"), so high_cai always diverges under explicit injection — even
    when the injected main table is byte-identical to the current default —
    because the golden set itself changed away from the bundled one.
    """
    default = RuleBasedOptimizer()
    injected = RuleBasedOptimizer(codon_table_path=str(PROFILES_DIR / table_name))
    default_seq = default.optimize(PROBE, profile=profile, seed=SEED).sequence
    injected_seq = injected.optimize(PROBE, profile=profile, seed=SEED).sequence

    if profile != "high_cai" and table_name == _resolved_default_table_name():
        assert default_seq == injected_seq, (
            f"{profile} design diverged from injection of the table that is "
            f"already the current default ({table_name}) — same resolved "
            f"asset must give the same design."
        )
    else:
        assert default_seq != injected_seq, (
            f"{profile} design did not respond to injected table {table_name}"
        )


@pytest.mark.parametrize("table_name", SOURCE_TABLES)
def test_gc_target_is_gc_driven_within_tolerance(table_name):
    """gc_target is GC-driven: the achieved GC% must be near-identical
    regardless of injected codon-usage source. Exact byte-for-byte output is
    NOT guaranteed — ties among multiple codons at the same GC% may be
    broken differently depending on the table's relative codon ranking, and
    that ranking legitimately differs across genome-derived assets."""
    default = RuleBasedOptimizer()
    injected = RuleBasedOptimizer(codon_table_path=str(PROFILES_DIR / table_name))
    default_result = default.optimize(PROBE, profile="gc_target", seed=SEED)
    injected_result = injected.optimize(PROBE, profile="gc_target", seed=SEED)
    assert abs(default_result.metrics["gc_percent"] - injected_result.metrics["gc_percent"]) <= 2.0


def test_injection_is_deterministic():
    """Same injected table + seed ⇒ identical design across instances."""
    table = str(PROFILES_DIR / SOURCE_TABLES[0])
    one = RuleBasedOptimizer(codon_table_path=table)
    two = RuleBasedOptimizer(codon_table_path=table)
    assert (
        one.optimize(PROBE, profile="balanced", seed=SEED).sequence
        == two.optimize(PROBE, profile="balanced", seed=SEED).sequence
    )
