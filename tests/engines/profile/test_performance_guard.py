"""Performance guard tests for profile engine."""

from __future__ import annotations

import statistics
import sys
import time
from pathlib import Path

# Add project src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from factorforge.engines.profile.optimizer import RuleBasedOptimizer


def _load_sample_protein() -> str:
    path = Path(__file__).resolve().parents[3] / "examples" / "sample_protein.fasta"
    return "".join(
        line.strip() for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith(">")
    )


def test_fast_scan_mode_faster_than_full_on_long_input() -> None:
    """Fast scan mode should outperform full scanning on long inputs."""
    optimizer = RuleBasedOptimizer()
    base = _load_sample_protein()
    long_protein = (base * ((4000 // len(base)) + 1))[:4000]

    # Warmup to reduce one-time import/cache effects.
    optimizer.optimize(long_protein, profile="balanced", scan_mode="full")
    optimizer.optimize(long_protein, profile="balanced", scan_mode="fast")

    full_times = []
    fast_times = []
    for _ in range(3):
        t0 = time.perf_counter()
        optimizer.optimize(long_protein, profile="balanced", scan_mode="full")
        full_times.append(time.perf_counter() - t0)

        t1 = time.perf_counter()
        optimizer.optimize(long_protein, profile="balanced", scan_mode="fast")
        fast_times.append(time.perf_counter() - t1)

    assert statistics.median(fast_times) < statistics.median(full_times)


def test_balanced_profile_p95_under_2s() -> None:
    """Guard against major regressions in hot optimization path."""
    optimizer = RuleBasedOptimizer()
    protein = _load_sample_protein()

    durations = []
    for _ in range(20):
        t0 = time.perf_counter()
        optimizer.optimize(protein, profile="balanced", scan_mode="full")
        durations.append(time.perf_counter() - t0)

    durations_sorted = sorted(durations)
    p95_idx = max(0, int(len(durations_sorted) * 0.95) - 1)
    p95 = durations_sorted[p95_idx]
    assert p95 < 2.0
