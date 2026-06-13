# benchmarks/ablation/run_ablation.py
"""Main ablation runner.

Usage:
    python -m benchmarks.ablation.run_ablation [--limit N] [--skip-run]

--limit N   : cap sequences (for quick validation)
--skip-run  : skip L2/L3/L4 new runs (extract+summary only, requires existing ablation_results.csv)
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
import pandas as pd
import yaml

from benchmarks.config import load_benchmark_config
from benchmarks.scoring import score_cds, canonical_multi_constraint_pass

ROOT = Path(__file__).resolve().parents[2]
SPEC_PATH = ROOT / "benchmarks" / "ablation" / "ablation_spec.yaml"
OUT_DIR = ROOT / "benchmarks" / "results" / "v3.2.0" / "ablation"
OUT_CSV = OUT_DIR / "ablation_results.csv"
OUT_SUMMARY = OUT_DIR / "ablation_summary.json"


def _read_fasta(p: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def extract_existing_layers(
    formal_df: pd.DataFrame,
    ablation_layer: str,
    method_name: str,
    enabled_constraints: dict | None = None,
) -> pd.DataFrame:
    """Extract rows for a given method from the formal benchmark CSV."""
    df = formal_df[formal_df["method"] == method_name].copy()
    df["ablation_layer"] = ablation_layer
    df["ablation_source"] = "existing_csv"
    df["enabled_constraints_json"] = json.dumps(enabled_constraints or {})
    return df


def _load_condition_fn(module_path: str, fn_name: str):
    import importlib
    mod = importlib.import_module(module_path)
    return getattr(mod, fn_name)


def run_new_layer(
    proteins: dict[str, str],
    ablation_layer: str,
    layer_spec: dict,
    cfg,
    progress_every: int = 500,
) -> pd.DataFrame:
    """Run a new ablation condition on all proteins."""
    fn = _load_condition_fn(layer_spec["condition_module"], layer_spec["condition_fn"])
    max_attempts = layer_spec.get("max_attempts", 50)
    method_name = layer_spec["name"]
    enabled_constraints_str = json.dumps(layer_spec.get("enabled_constraints", {}))
    rows = []
    items = list(proteins.items())
    total = len(items)
    print(f"[ablation] {ablation_layer} ({method_name}): running {total} sequences", flush=True)
    for i, (sid, protein) in enumerate(items):
        if (i + 1) % progress_every == 0 or (i + 1) == total:
            print(f"[ablation] {ablation_layer} {i + 1}/{total}", flush=True)
        t0 = time.perf_counter()
        try:
            if "max_attempts" in layer_spec:
                cds = fn(protein, seed=320, max_attempts=max_attempts)
            else:
                cds = fn(protein)
            rt = time.perf_counter() - t0
            row = score_cds(method_name, "ablation", sid, protein, cds, cfg, rt)
            row["status"] = "ok"
            row["error_type"] = None
            row["error_message"] = None
        except Exception as exc:
            rt = time.perf_counter() - t0
            row = {
                "method": method_name, "method_type": "ablation",
                "sequence_id": sid,
                "aa_identity": None, "internal_stop_count": None,
                "invalid_codon_count": None, "length_multiple_of_three": None,
                "cai": None, "gc_percent": None, "gc_in_target_range": None,
                "forbidden_type_iis_site_count": None,
                "biological_pass": False, "assembly_pass": False,
                "multi_constraint_pass": False, "runtime_seconds": round(rt, 6),
                "status": "error", "error_type": type(exc).__name__,
                "error_message": str(exc)[:300],
            }
        row["replicate"] = 1
        row["seed"] = 320
        row["ablation_layer"] = ablation_layer
        row["ablation_source"] = "new_run"
        row["enabled_constraints_json"] = enabled_constraints_str
        rows.append(row)
    return pd.DataFrame(rows)


def build_ablation_summary(
    df: pd.DataFrame,
    spec: dict,
    formal_run_id: str,
    formal_summary_sha256: str,
    spec_sha256: str,
    input_fasta_sha256: str,
) -> dict:
    """Build ablation_summary.json from the merged ablation DataFrame."""
    layers_summary = {}
    for layer_name in ["L0", "L1", "L2", "L3", "L4", "L5"]:
        sub = df[df["ablation_layer"] == layer_name]
        if sub.empty:
            continue
        ok = sub[sub["status"] == "ok"] if "status" in sub.columns else sub
        if layer_name == "L0":
            ok = ok.groupby("sequence_id")[
                ["multi_constraint_pass", "biological_pass", "assembly_pass",
                 "gc_in_target_range", "cai", "gc_percent",
                 "forbidden_type_iis_site_count"]
            ].mean().reset_index()
        n_seqs = len(ok["sequence_id"].unique()) if "sequence_id" in ok.columns else len(ok)
        layers_summary[layer_name] = {
            "method_name": spec["layers"][layer_name]["name"],
            "source": spec["layers"][layer_name]["source"],
            "enabled_constraints": spec["layers"][layer_name]["enabled_constraints"],
            "n_sequences": int(n_seqs),
            "multi_constraint_pass_rate": round(float(ok["multi_constraint_pass"].mean()), 4),
            "biological_pass_rate": round(float(ok["biological_pass"].mean()), 4),
            "assembly_pass_rate": round(float(ok["assembly_pass"].mean()), 4),
            "gc_in_range_rate": round(float(ok["gc_in_target_range"].mean()), 4),
            "mean_cai": round(float(ok["cai"].mean()), 4),
        }
    return {
        "analysis_type": "constraint_ablation",
        "scoring_contract_version": "v1.1",
        "multi_constraint_pass_definition": "biological_pass AND assembly_pass AND gc_in_target_range",
        "source_formal_run_id": formal_run_id,
        "source_formal_summary_sha256": formal_summary_sha256,
        "ablation_spec_sha256": spec_sha256,
        "input_fasta_sha256": input_fasta_sha256,
        "dataset_n": spec["ablation"]["dataset_n"],
        "seed": spec["ablation"]["random_seed"],
        "layers": layers_summary,
    }


def main(limit: int | None = None, skip_run: bool = False) -> None:
    spec = yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))
    cfg = load_benchmark_config()

    formal_csv = ROOT / spec["ablation"]["source_formal_results_csv"]
    proteins_fasta = ROOT / spec["ablation"]["input_fasta"]
    frozen_summary = ROOT / spec["ablation"]["source_formal_summary"]

    formal_df = pd.read_csv(formal_csv)
    proteins = _read_fasta(proteins_fasta)
    if limit:
        proteins = dict(list(proteins.items())[:limit])

    # Filter formal_df to same sequence_id subset as proteins (critical for --limit)
    sequence_ids = set(proteins.keys())
    formal_df = formal_df[formal_df["sequence_id"].isin(sequence_ids)].copy()

    l0_constraints = spec["layers"]["L0"]["enabled_constraints"]
    l1_constraints = spec["layers"]["L1"]["enabled_constraints"]
    l5_constraints = spec["layers"]["L5"]["enabled_constraints"]
    l0 = extract_existing_layers(formal_df, "L0", "random_synonymous", l0_constraints)
    l1 = extract_existing_layers(formal_df, "L1", "greedy_cai", l1_constraints)
    l5 = extract_existing_layers(formal_df, "L5", "factorforge_assembly_friendly", l5_constraints)

    if skip_run:
        existing = pd.read_csv(OUT_CSV)
        l2 = existing[existing["ablation_layer"] == "L2"]
        l3 = existing[existing["ablation_layer"] == "L3"]
        l4 = existing[existing["ablation_layer"] == "L4"]
    else:
        l2 = run_new_layer(proteins, "L2", spec["layers"]["L2"], cfg)
        l3 = run_new_layer(proteins, "L3", spec["layers"]["L3"], cfg)
        l4 = run_new_layer(proteins, "L4", spec["layers"]["L4"], cfg)

    merged = pd.concat([l0, l1, l2, l3, l4, l5], ignore_index=True)
    # Recompute multi_constraint_pass from primitive columns (scoring_contract v1.1).
    # This corrects stale values in CSVs produced before the gc_in_target_range fix.
    merged["multi_constraint_pass"] = canonical_multi_constraint_pass(
        merged, gc_min=cfg.gc_min, gc_max=cfg.gc_max
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_CSV, index=False)
    print(f"[ablation] wrote {len(merged)} rows → {OUT_CSV}", flush=True)

    spec_sha256 = hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest()
    input_sha256 = hashlib.sha256(proteins_fasta.read_bytes()).hexdigest() if proteins_fasta.exists() else ""
    if frozen_summary.exists():
        frozen_data = json.loads(frozen_summary.read_text(encoding="utf-8"))
        formal_run_id = frozen_data.get("run_id", "")
        formal_sha256 = hashlib.sha256(frozen_summary.read_bytes()).hexdigest()
    else:
        formal_run_id, formal_sha256 = "", ""
    try:
        _git = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        _git = None

    summary = build_ablation_summary(merged, spec, formal_run_id, formal_sha256, spec_sha256, input_sha256)
    summary["git_commit"] = _git
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[ablation] wrote summary → {OUT_SUMMARY}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--skip-run", action="store_true")
    args = ap.parse_args()
    main(limit=args.limit, skip_run=args.skip_run)
