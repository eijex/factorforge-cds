"""Compare Run 4B smoke result directories against acceptance criteria."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _bool_pass_rate(rows: list[dict[str, str]], field: str) -> float:
    if not rows:
        return 0.0
    return sum(1 for row in rows if str(row[field]).lower() == "true") / len(rows)


def _load_result_dir(result_dir: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    comparison_path = result_dir / "run4_comparison.csv"
    loss_path = result_dir / "run4_loss_log.csv"
    summary_path = result_dir / "run4_summary.json"
    if not comparison_path.exists():
        raise FileNotFoundError(f"Missing comparison CSV: {comparison_path}")
    if not loss_path.exists():
        raise FileNotFoundError(f"Missing loss log CSV: {loss_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    return _read_csv(comparison_path), _read_csv(loss_path), summary


def summarize_result(result_dir: str | Path) -> dict[str, Any]:
    """Summarize one Run 4B smoke result directory."""
    path = Path(result_dir)
    comparison, loss, raw_summary = _load_result_dir(path)
    if not comparison:
        raise ValueError(f"No comparison rows in {path}")
    if not loss:
        raise ValueError(f"No loss rows in {path}")

    last_step = max(int(row["step"]) for row in loss)
    last_eval = [row for row in loss if int(row["step"]) == last_step]
    v2_cai = [float(row["v2_cai"]) for row in comparison]
    v3_cai = [float(row["v3_cai"]) for row in comparison]
    ratio = [v3 / v2 if v2 > 0 else 1.0 for v2, v3 in zip(v2_cai, v3_cai)]
    validator_pass_rate = _bool_pass_rate(comparison, "v3_validator_pass")
    identity_values = [float(row["v3_amino_acid_identity"]) for row in comparison]
    v2_gc = [float(row["v2_gc_global"]) for row in comparison]
    v3_gc = [float(row["v3_gc_global"]) for row in comparison]

    summary = {
        "result_dir": str(path),
        "decoded_count": int(len(comparison)),
        "validator_pass_rate": validator_pass_rate,
        "amino_acid_identity_min": min(identity_values),
        "v2_cai_mean": _mean(v2_cai),
        "v3_cai_mean": _mean(v3_cai),
        "v3_to_v2_cai_ratio_min": min(ratio),
        "v3_to_v2_cai_ratio_mean": _mean(ratio),
        "v2_global_gc_mean": _mean(v2_gc),
        "v3_global_gc_mean": _mean(v3_gc),
        "v3_global_gc_min": min(v3_gc),
        "v3_global_gc_max": max(v3_gc),
        "expected_gc_mean": _mean([float(row["expected_GC"]) for row in last_eval]),
        "bounded_gc_penalty_mean": _mean(
            [float(row["bounded_GC_penalty"]) for row in last_eval]
        ),
        "raw_summary": raw_summary,
    }
    criteria = {
        "validator_pass_rate_100pct": summary["validator_pass_rate"] == 1.0,
        "amino_acid_identity_100pct": summary["amino_acid_identity_min"] == 1.0,
        "decoded_global_gc_mean_ge_40pct": summary["v3_global_gc_mean"] >= 40.0,
        "decoded_global_gc_max_le_60pct": summary["v3_global_gc_max"] <= 60.0,
        "expected_gc_ge_0_40": summary["expected_gc_mean"] >= 0.40,
        "v3_v2_cai_ratio_mean_ge_0_97": summary["v3_to_v2_cai_ratio_mean"] >= 0.97,
    }
    summary["acceptance_criteria"] = criteria
    summary["verdict"] = "pass" if all(criteria.values()) else "fail"
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare Run 4B smoke result directories.")
    parser.add_argument("result_dirs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    summaries = [summarize_result(path) for path in args.result_dirs]
    payload = {"results": summaries}
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
