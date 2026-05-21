"""Run reproducible v2 baseline metrics for v3-alpha benchmark proteins."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TESTS = ROOT / "tests"
for path in (SRC, TESTS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from factorforge.engines.v3.inference.v2_adapter import optimize_with_v2
from fixtures.benchmark_proteins import BENCHMARK_PROTEINS


OUTPUT_DIR = ROOT / "experiments" / "results" / "metrics"


def _row(name: str, result: dict[str, Any]) -> dict[str, Any]:
    validator = result["validator"]
    return {
        "name": name,
        "engine": result["engine"],
        "profile": result["metadata"]["profile"],
        "protein_length": len(result["protein_sequence"]),
        "dna_length": len(result["dna_sequence"]),
        "cai": result["metrics"]["cai"],
        "gc": result["metrics"]["gc"],
        "amino_acid_identity": validator["amino_acid_identity"],
        "internal_stop_count": validator["internal_stop_count"],
        "invalid_codon_count": validator["invalid_codon_count"],
        "validator_passed": validator["passed"],
    }


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []

    for name, protein in BENCHMARK_PROTEINS.items():
        result = optimize_with_v2(protein, options={"profile": "high_cai", "scan_mode": "fast"})
        results.append({"name": name, **result})
        rows.append(_row(name, result))

    csv_path = OUTPUT_DIR / "v2_baseline_metrics.csv"
    json_path = OUTPUT_DIR / "v2_baseline_metrics.json"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    json_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    for row in rows:
        print(
            f"{row['name']}: CAI={row['cai']:.3f}, GC={row['gc']:.1f}%, "
            f"AA identity={row['amino_acid_identity']:.3f}, passed={row['validator_passed']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

