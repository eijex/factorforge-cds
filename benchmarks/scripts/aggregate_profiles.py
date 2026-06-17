"""Aggregate A/B/C1/C2 frozen benchmark summaries into comparison tables."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
for _path in (ROOT, SRC):
    _path_str = str(_path)
    if _path_str not in sys.path:
        sys.path.insert(0, _path_str)

from factorforge.analysis.metrics import load_codon_usage_table  # noqa: E402


OUT_DIR = ROOT / "reproducibility" / "benchmark_v3.2.2"
PROFILES = {
    "A": {
        "profile_id": "legacy_packaged",
        "summary": "benchmark_summary_profile_legacy_packaged.frozen.json",
        "codon_table": ROOT / "src" / "factorforge" / "data" / "nbenthamiana_codons.json",
    },
    "B": {
        "profile_id": "qld183_v103_derived",
        "summary": "benchmark_summary_profile_qld183_v103_derived.frozen.json",
        "codon_table": ROOT / "src" / "factorforge" / "data" / "profiles" / "qld183_v103_derived_codons.json",
    },
    "C1": {
        "profile_id": "nbev11_cds_all_derived",
        "summary": "benchmark_summary_profile_nbev11_cds_all_derived.frozen.json",
        "codon_table": ROOT / "src" / "factorforge" / "data" / "profiles" / "nbev11_cds_all_derived_codons.json",
    },
    "C2": {
        "profile_id": "nbev11_cds_hc_derived",
        "summary": "benchmark_summary_profile_nbev11_cds_hc_derived.frozen.json",
        "codon_table": ROOT / "src" / "factorforge" / "data" / "profiles" / "nbev11_cds_hc_derived_codons.json",
    },
}


def _load_summaries() -> dict[str, dict] | None:
    missing = [OUT_DIR / item["summary"] for item in PROFILES.values() if not (OUT_DIR / item["summary"]).exists()]
    if missing:
        print("Missing frozen profile summaries:")
        for path in missing:
            print(f"- {path}")
        return None
    return {
        label: json.loads((OUT_DIR / item["summary"]).read_text(encoding="utf-8"))
        for label, item in PROFILES.items()
    }


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        return float("nan")
    mean_left = sum(left) / len(left)
    mean_right = sum(right) / len(right)
    numerator = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right))
    denom_left = math.sqrt(sum((a - mean_left) ** 2 for a in left))
    denom_right = math.sqrt(sum((b - mean_right) ** 2 for b in right))
    if denom_left == 0 or denom_right == 0:
        return float("nan")
    return numerator / (denom_left * denom_right)


def _load_weight_tables() -> dict[str, dict[str, float]]:
    return {
        label: load_codon_usage_table(item["codon_table"]).codon_weights
        for label, item in PROFILES.items()
    }


def _combined_rows(summaries: dict[str, dict]) -> list[dict]:
    rows = []
    for label, summary in summaries.items():
        profile_id = PROFILES[label]["profile_id"]
        for method, metrics in summary["methods"].items():
            for metric, value in metrics.items():
                if isinstance(value, int | float):
                    rows.append(
                        {
                            "profile_id": profile_id,
                            "method": method,
                            "metric": metric,
                            "value": value,
                        }
                    )
    return rows


def _multi_constraint_rows(summaries: dict[str, dict]) -> list[dict]:
    methods = sorted({method for summary in summaries.values() for method in summary["methods"]})
    rows = []
    for method in methods:
        row = {"method": method}
        for label, summary in summaries.items():
            row[label] = summary["methods"].get(method, {}).get("pass_rate_multi_constraint")
        rows.append(row)
    return rows


def _cai_rows(summaries: dict[str, dict]) -> list[dict]:
    methods = sorted({method for summary in summaries.values() for method in summary["methods"]})
    rows = []
    for method in methods:
        values = {
            label: summary["methods"].get(method, {}).get("mean_cai")
            for label, summary in summaries.items()
        }
        baseline = values["A"]
        rows.append(
            {
                "method": method,
                "mean_cai_A": values["A"],
                "mean_cai_B": values["B"],
                "mean_cai_C1": values["C1"],
                "mean_cai_C2": values["C2"],
                "shift_B_vs_A": None if baseline is None or values["B"] is None else round(values["B"] - baseline, 4),
                "shift_C1_vs_A": None if baseline is None or values["C1"] is None else round(values["C1"] - baseline, 4),
                "shift_C2_vs_A": None if baseline is None or values["C2"] is None else round(values["C2"] - baseline, 4),
            }
        )
    return rows


def _correlation_rows(weights: dict[str, dict[str, float]]) -> list[dict]:
    codon_order = list(weights["A"].keys())
    rows = []
    for left_label in PROFILES:
        row = {"profile": left_label}
        left_values = [weights[left_label][codon] for codon in codon_order]
        for right_label in PROFILES:
            right_values = [weights[right_label][codon] for codon in codon_order]
            row[right_label] = round(_pearson(left_values, right_values), 6)
        rows.append(row)
    return rows


def _rank_map(weights: dict[str, float]) -> dict[str, int]:
    aa_to_codons: dict[str, list[str]] = {}
    table = load_codon_usage_table(PROFILES["A"]["codon_table"])
    for codon, aa in table.codon_to_aa.items():
        if aa != "*" and codon in weights:
            aa_to_codons.setdefault(aa, []).append(codon)

    ranks = {}
    for codons in aa_to_codons.values():
        ranked = sorted(codons, key=lambda codon: (-weights[codon], codon))
        for rank, codon in enumerate(ranked, start=1):
            ranks[codon] = rank
    return ranks


def _preferred_rank_rows(weights: dict[str, dict[str, float]]) -> list[dict]:
    legacy_table = load_codon_usage_table(PROFILES["A"]["codon_table"])
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in legacy_table.codon_to_aa.items():
        if aa != "*" and codon in legacy_table.codon_weights:
            aa_to_codons.setdefault(aa, []).append(codon)

    rank_maps = {label: _rank_map(table_weights) for label, table_weights in weights.items()}
    rows = []
    for aa in sorted(aa_to_codons):
        codons = aa_to_codons[aa]
        if len(codons) < 2:
            continue
        for codon in codons:
            rows.append(
                {
                    "aa": aa,
                    "codon": codon,
                    "rank_A": rank_maps["A"][codon],
                    "rank_B": rank_maps["B"][codon],
                    "rank_C1": rank_maps["C1"][codon],
                    "rank_C2": rank_maps["C2"][codon],
                    "rank_shift_C1_vs_A": rank_maps["C1"][codon] - rank_maps["A"][codon],
                }
            )
    return rows


def main() -> int:
    summaries = _load_summaries()
    if summaries is None:
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    combined = _combined_rows(summaries)
    (OUT_DIR / "combined_profile_comparison.json").write_text(
        json.dumps(combined, indent=2) + "\n",
        encoding="utf-8",
    )
    _write_tsv(
        OUT_DIR / "combined_profile_comparison.tsv",
        combined,
        ["profile_id", "method", "metric", "value"],
    )

    multi_rows = _multi_constraint_rows(summaries)
    _write_tsv(OUT_DIR / "multi_constraint_pass_by_profile.tsv", multi_rows, ["method", "A", "B", "C1", "C2"])

    cai_rows = _cai_rows(summaries)
    _write_tsv(
        OUT_DIR / "cai_shift_by_profile.tsv",
        cai_rows,
        [
            "method",
            "mean_cai_A",
            "mean_cai_B",
            "mean_cai_C1",
            "mean_cai_C2",
            "shift_B_vs_A",
            "shift_C1_vs_A",
            "shift_C2_vs_A",
        ],
    )

    weights = _load_weight_tables()
    _write_tsv(
        OUT_DIR / "codon_weight_correlation_matrix.tsv",
        _correlation_rows(weights),
        ["profile", "A", "B", "C1", "C2"],
    )
    _write_tsv(
        OUT_DIR / "preferred_codon_rank_shift.tsv",
        _preferred_rank_rows(weights),
        ["aa", "codon", "rank_A", "rank_B", "rank_C1", "rank_C2", "rank_shift_C1_vs_A"],
    )

    print("| profile | dataset_n | source_profile_id |")
    print("|---|---:|---|")
    for label, summary in summaries.items():
        print(f"| {label} | {summary.get('dataset_n')} | {summary.get('source_profile_id')} |")
    print(f"Aggregated outputs written to: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
