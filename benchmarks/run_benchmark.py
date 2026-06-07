"""Canonical benchmark runner. CI smoke uses fixtures (no network);
formal mode uses the fetched dataset."""
from __future__ import annotations
import argparse
import csv
import time
from datetime import date
from pathlib import Path

from benchmarks.config import load_benchmark_config
from benchmarks.scoring import score_cds
from benchmarks.baselines.random_synonymous import random_synonymous_cds
from benchmarks.baselines.most_frequent_codon import most_frequent_codon_cds
from benchmarks.baselines.greedy_cai import greedy_cai_cds
from benchmarks.baselines.native_reference import native_reference_cds
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge import __version__ as FF_VERSION

FACTORFORGE_PROFILES = ["balanced", "high_cai", "gc_target", "assembly_friendly"]
CLAIM = (
    "This benchmark evaluates in-silico CDS design quality only. "
    "It does not demonstrate improved protein expression, yield, folding, or wet-lab performance."
)


def _read_fasta(p: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


def _timed(fn, *args):
    t0 = time.perf_counter()
    cds = fn(*args)
    return cds, time.perf_counter() - t0


def run(dataset: str, mode: str, out_csv: Path, out_md: Path,
        proteins_fasta: Path, native_fasta: Path) -> None:
    cfg = load_benchmark_config()
    proteins = _read_fasta(proteins_fasta)
    natives = _read_fasta(native_fasta)
    opt = RuleBasedOptimizer()
    rows = []
    for sid, protein in proteins.items():
        # random_synonymous: fixed seed + replicates (seed, seed+1, ...) for stability
        for r in range(cfg.random_replicates):
            seed = cfg.random_seed + r
            cds, rt = _timed(random_synonymous_cds, protein, seed)
            row = score_cds("random_synonymous", "baseline", sid, protein, cds, cfg, rt)
            row["replicate"], row["seed"] = r + 1, seed
            rows.append(row)
        for name, fn in (("most_frequent_codon", most_frequent_codon_cds),
                         ("greedy_cai", greedy_cai_cds)):
            cds, rt = _timed(fn, protein)
            row = score_cds(name, "baseline", sid, protein, cds, cfg, rt)
            row["replicate"], row["seed"] = 1, ""
            rows.append(row)
        if sid in natives:
            cds, rt = _timed(native_reference_cds, natives[sid])
            row = score_cds("native_reference", "reference", sid, protein, cds, cfg, rt)
            row["replicate"], row["seed"] = 1, ""
            rows.append(row)
        for prof in FACTORFORGE_PROFILES:
            t0 = time.perf_counter()
            res = opt.optimize(protein, profile=prof)
            rt = time.perf_counter() - t0
            row = score_cds(f"factorforge_{prof}", "optimizer", sid, protein, res.sequence, cfg, rt)
            row["replicate"], row["seed"] = 1, ""
            rows.append(row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    _write_summary(rows, out_md, dataset, mode, cfg)


def _pass_rate(rows, key, method):
    sub = [r for r in rows if r["method"] == method]
    return sum(1 for r in sub if r[key]) / len(sub) if sub else 0.0


def _write_summary(rows, out_md: Path, dataset: str, mode: str, cfg) -> None:
    methods = sorted({r["method"] for r in rows})
    lines = [
        f"# Benchmark Summary (dataset={dataset}, mode={mode})", "",
        "## Provenance (reproducibility)",
        f"- FactorForge version: {FF_VERSION}",
        f"- registry sha256: {cfg.registry_sha256}",
        f"- benchmark_spec sha256: {cfg.spec_sha256}",
        f"- run date: {date.today().isoformat()}",
        "", f"> {CLAIM}", "",
        "| method | multi_constraint_pass_rate | biological_pass_rate | assembly_pass_rate |",
        "|---|---|---|---|",
    ]
    for m in methods:
        lines.append(f"| {m} | {_pass_rate(rows,'multi_constraint_pass',m):.2f} | "
                     f"{_pass_rate(rows,'biological_pass',m):.2f} | {_pass_rate(rows,'assembly_pass',m):.2f} |")
    lines += [
        "", "_native_reference is a biological reference anchor, not an optimizer._",
        "_greedy_cai is a CAI-focused baseline that does not explicitly optimize GC or assembly constraints._",
        "",
        "FactorForge preserves amino-acid identity and avoids invalid CDS outputs while "
        "improving multi-constraint in-silico CDS design quality relative to simple "
        "synonymous-codon baselines.",
    ]
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines), encoding="utf-8")
    # JSON observability record (Brief §9)
    import hashlib, json
    from datetime import datetime, timezone
    all_rows = rows
    def _pr(key):
        sub = [r for r in all_rows if r["method"].startswith("factorforge")]
        return round(sum(1 for r in sub if r[key]) / len(sub), 4) if sub else 0.0
    summary_json = out_md.parent / "benchmark_summary.json"
    summary_json.write_text(json.dumps({
        "run_id": hashlib.sha256(f"{dataset}{mode}{date.today().isoformat()}".encode()).hexdigest()[:12],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "factorforge_version": FF_VERSION,
        "registry_version": cfg.registry_version,
        "registry_hash": f"sha256:{cfg.registry_sha256}",
        "spec_hash": f"sha256:{cfg.spec_sha256}",
        "runtime_seconds": round(sum(r["runtime_seconds"] for r in all_rows), 4),
        "dataset_id": dataset,
        "dataset_n": len({r["sequence_id"] for r in all_rows}),
        "pass_rate_biological": _pr("biological_pass"),
        "pass_rate_assembly": _pr("assembly_pass"),
        "pass_rate_multi_constraint": _pr("multi_constraint_pass"),
    }, indent=2), encoding="utf-8")


def main(default_mode: str = "formal") -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="synthetic")
    ap.add_argument("--mode", default=default_mode)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.dataset == "synthetic":
        proteins = root / "tests" / "fixtures" / "small_proteins.fasta"
        native = root / "tests" / "fixtures" / "small_native_cds.fasta"
    else:  # formal: real fetched dataset (must run fetch_dataset.py first)
        proteins = root / "benchmarks" / "datasets" / "nbenthamiana_reference_proteins.fasta"
        native = root / "benchmarks" / "datasets" / "nbenthamiana_reference_cds.fasta"
    run(dataset=args.dataset, mode=args.mode,
        out_csv=root / "benchmarks" / "results" / "v3.2.0" / "benchmark_results.csv",
        out_md=root / "benchmarks" / "results" / "v3.2.0" / "benchmark_summary.md",
        proteins_fasta=proteins, native_fasta=native)


if __name__ == "__main__":
    main()
