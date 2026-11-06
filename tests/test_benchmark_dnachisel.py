"""Unit tests for DNA Chisel vs FactorForge comparative benchmark."""

import pytest

pytest.importorskip("dnachisel", reason="Install the benchmark extra for cross-tool evaluation")

from benchmarks.benchmark_dnachisel_vs_factorforge import run_benchmark

ADALIMUMAB_VH = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
)


def test_dnachisel_vs_factorforge_benchmark_execution() -> None:
    targets = {"Adalimumab_VH": ADALIMUMAB_VH}

    records = run_benchmark(
        targets=targets,
        gc_min=0.40,
        gc_max=0.46,
        host_name="nbenthamiana",
        verbose=False,
    )

    assert len(records) == 2  # 1 FactorForge + 1 DNA Chisel
    ff_rec = next(r for r in records if r["engine"] == "FactorForge DP")
    dc_rec = next(r for r in records if r["engine"] == "DNA Chisel")

    # FactorForge assertions: exact DP guarantees constraint feasibility
    assert ff_rec["status"] == "FEASIBLE"
    assert ff_rec["motif_violations"] == 0
    assert 40.0 <= ff_rec["gc_percent"] <= 46.0
    assert ff_rec["optimality"] == "GLOBAL_EXACT"
    assert ff_rec["cai"] >= 0.88

    # DNA Chisel assertions: stochastic heuristic feasibility
    assert dc_rec["status"] == "FEASIBLE"
    assert dc_rec["motif_violations"] == 0
    assert 40.0 <= dc_rec["gc_percent"] <= 46.0
    assert dc_rec["cai"] >= 0.88

    # Algorithmic boundary: FactorForge incorporates a 5' AU initiation ramp (codons 1-15)
    # in its discrete objective function, so raw unpenalized CAI is comparable within tolerance.
    assert ff_rec["cai"] >= dc_rec["cai"] - 0.01, (
        f"FactorForge CAI ({ff_rec['cai']}) should be comparable to or outperform DNA Chisel CAI ({dc_rec['cai']})"
    )
