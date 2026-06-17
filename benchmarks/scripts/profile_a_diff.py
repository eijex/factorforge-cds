"""P8: Compare fresh Profile A benchmark run against frozen v0.5.1 summary."""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "reproducibility" / "benchmark_v0.5.1" / "data" / "benchmark_summary.frozen.json"
FRESH  = ROOT / "benchmarks" / "results" / "v3.2.0" / "benchmark_summary.json"
OUT_TSV  = ROOT / "reproducibility" / "benchmark_v0.5.1" / "data" / "profile_A_reproduction_diff.tsv"
OUT_JSON = ROOT / "reproducibility" / "benchmark_v0.5.1" / "data" / "profile_A_reproduction_diff.json"

METRICS = [
    "pass_rate_multi_constraint",
    "pass_rate_biological",
    "pass_rate_assembly",
    "gc_in_range_rate",
    "mean_cai",
]

frozen = json.loads(FROZEN.read_text(encoding="utf-8"))
fresh  = json.loads(FRESH.read_text(encoding="utf-8"))

rows = []
for method, fz in frozen["methods"].items():
    fr = fresh["methods"].get(method, {})
    for metric in METRICS:
        fz_val = fz.get(metric)
        fr_val = fr.get(metric)
        diff = round(fr_val - fz_val, 6) if (fz_val is not None and fr_val is not None) else None
        match = (diff == 0.0) if diff is not None else None
        rows.append({
            "method": method,
            "metric": metric,
            "frozen": fz_val,
            "fresh": fr_val,
            "diff": diff,
            "match": match,
        })

verdict = "PASS" if all(r["match"] for r in rows if r["match"] is not None) else "DIVERGED"
diverged = [r for r in rows if r["match"] is False]

OUT_TSV.write_text(
    "\t".join(["method", "metric", "frozen", "fresh", "diff", "match"]) + "\n" +
    "\n".join("\t".join(str(r[k]) for k in ["method", "metric", "frozen", "fresh", "diff", "match"]) for r in rows),
    encoding="utf-8",
)
OUT_JSON.write_text(
    json.dumps({"verdict": verdict, "diverged_count": len(diverged), "rows": rows}, indent=2),
    encoding="utf-8",
)

print(f"P8 verdict: {verdict}")
print(f"Diverged rows: {len(diverged)}")
if diverged:
    for r in diverged:
        print(f"  {r['method']} / {r['metric']}: frozen={r['frozen']} fresh={r['fresh']} diff={r['diff']}")
print(f"TSV: {OUT_TSV}")
print(f"JSON: {OUT_JSON}")
