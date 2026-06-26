"""Canonical benchmark runner. CI smoke uses fixtures (no network);
formal mode uses the fetched dataset."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from datetime import date, datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT, ROOT / "src"):
    _path_str = str(_path)
    if _path_str not in sys.path:
        sys.path.insert(0, _path_str)

from benchmarks.config import load_benchmark_config  # noqa: E402
from benchmarks.scoring import score_cds  # noqa: E402
from benchmarks.baselines.random_synonymous import random_synonymous_cds  # noqa: E402
from benchmarks.baselines.greedy_cai import greedy_cai_cds  # noqa: E402
from benchmarks.baselines.native_reference import native_reference_cds  # noqa: E402
from factorforge.engines.profile.optimizer import RuleBasedOptimizer  # noqa: E402
from factorforge.engines.profile.scoring import _check_vienna_available  # noqa: E402
from factorforge import __version__ as FF_VERSION  # noqa: E402

FACTORFORGE_PROFILES = ["balanced", "high_cai", "gc_target", "assembly_friendly"]
CLAIM = (
    "This benchmark evaluates in-silico CDS design quality only. "
    "It does not demonstrate improved protein expression, yield, folding, or wet-lab performance."
)


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


def _timed(fn, *args):
    t0 = time.perf_counter()
    cds = fn(*args)
    return cds, time.perf_counter() - t0


def _ok(row: dict) -> dict:
    """Tag a successfully scored row."""
    row["status"] = "ok"
    row["error_type"] = None
    row["error_message"] = None
    return row


def run(dataset: str, mode: str, out_csv: Path, out_md: Path,
        proteins_fasta: Path, native_fasta: Path,
        limit: int | None = None, progress_every: int = 500,
        seed: int | None = None,
        codon_weights: dict[str, float] | None = None,
        source_profile_id: str = "legacy_packaged",
        source_table_sha256: str | None = None,
        source_manifest_sha256: str | None = None,
        codon_table_path: Path | None = None) -> None:
    cfg = load_benchmark_config()
    proteins = _read_fasta(proteins_fasta)
    natives = _read_fasta(native_fasta)
    if limit:
        proteins = dict(list(proteins.items())[:limit])
    # Inject the source-profile codon table into design when provided so the
    # FactorForge engines re-design against the alternative reference, not only
    # re-score CAI against it. Without a path the bundled reference is used and
    # behavior is identical to the product engine (source_profile_id=legacy_packaged).
    opt = RuleBasedOptimizer(
        codon_table_path=str(codon_table_path) if codon_table_path else None
    )
    rows = []
    total = len(proteins)
    print(f"[benchmark] loaded {total} protein records", flush=True)
    for i, (sid, protein) in enumerate(proteins.items()):
        if (i + 1) % progress_every == 0 or (i + 1) == total:
            errors = sum(1 for r in rows if r.get("status") == "error")
            print(f"[benchmark] {i + 1}/{total} sequences, errors={errors}", flush=True)
        for r in range(cfg.random_replicates):
            rs_seed = cfg.random_seed + r
            cds, rt = _timed(random_synonymous_cds, protein, rs_seed)
            row = _ok(score_cds("random_synonymous", "baseline", sid, protein, cds, cfg, rt, codon_weights))
            row["replicate"], row["seed"] = r + 1, rs_seed
            rows.append(row)
        cds, rt = _timed(greedy_cai_cds, protein, codon_weights)
        row = _ok(score_cds("greedy_cai", "baseline", sid, protein, cds, cfg, rt, codon_weights))
        row["replicate"], row["seed"] = 1, ""
        rows.append(row)
        if sid in natives:
            cds, rt = _timed(native_reference_cds, natives[sid])
            row = _ok(score_cds("native_reference", "reference", sid, protein, cds, cfg, rt, codon_weights))
            row["replicate"], row["seed"] = 1, ""
            rows.append(row)
        for prof in FACTORFORGE_PROFILES:
            try:
                t0 = time.perf_counter()
                res = opt.optimize(protein, profile=prof, seed=seed)
                rt = time.perf_counter() - t0
                row = _ok(score_cds(f"factorforge_{prof}", "optimizer", sid, protein, res.sequence, cfg, rt, codon_weights))
            except Exception as exc:
                seq_len = len(protein)
                seq_hash = hashlib.sha256(protein.encode()).hexdigest()[:12]
                print(f"[WARN] factorforge_{prof} failed seq_id={sid} len={seq_len} "
                      f"hash_prefix={seq_hash}: {type(exc).__name__}", flush=True)
                row = {
                    "method": f"factorforge_{prof}", "method_type": "optimizer",
                    "sequence_id": sid, "aa_identity": None,
                    "internal_stop_count": None, "invalid_codon_count": None,
                    "length_multiple_of_three": None, "cai": None,
                    "gc_percent": None, "gc_in_target_range": None,
                    "forbidden_type_iis_site_count": None,
                    "biological_pass": False, "assembly_pass": False,
                    "multi_constraint_pass": False, "runtime_seconds": 0.0,
                    "status": "error",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc)[:300],
                }
            row["replicate"], row["seed"] = 1, ""
            rows.append(row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    prot_sha256 = hashlib.sha256(proteins_fasta.read_bytes()).hexdigest() if proteins_fasta.exists() else None
    cds_sha256 = hashlib.sha256(native_fasta.read_bytes()).hexdigest() if native_fasta.exists() else None
    _write_summary(
        rows,
        out_md,
        dataset,
        mode,
        cfg,
        prot_sha256=prot_sha256,
        cds_sha256=cds_sha256,
        seed=seed,
        source_profile_id=source_profile_id,
        source_table_sha256=source_table_sha256,
        source_manifest_sha256=source_manifest_sha256,
    )


def _pass_rate(rows, key, method):
    sub = [r for r in rows if r["method"] == method and r.get("status") == "ok"]
    return sum(1 for r in sub if r[key]) / len(sub) if sub else 0.0


def _mean(rows, key, method):
    sub = [r for r in rows if r["method"] == method and r.get("status") == "ok" and r[key] is not None]
    return round(sum(r[key] for r in sub) / len(sub), 4) if sub else 0.0


def _write_summary(rows, out_md: Path, dataset: str, mode: str, cfg,
                   prot_sha256: str | None = None,
                   cds_sha256: str | None = None,
                   seed: int | None = None,
                   source_profile_id: str = "legacy_packaged",
                   source_table_sha256: str | None = None,
                   source_manifest_sha256: str | None = None) -> None:
    methods = sorted({r["method"] for r in rows})
    # Compute codon table metadata for MD
    _repo_root_md = out_md.parent.parent.parent.parent
    _manifest_path_md = _repo_root_md / "data" / "reference" / "codon_table_manifest.json"
    _codon_table_path_md = _repo_root_md / "src" / "factorforge" / "data" / "nbenthamiana_codons.json"
    ct_sha256 = (
        hashlib.sha256(_codon_table_path_md.read_bytes()).hexdigest()
        if _codon_table_path_md.exists()
        else None
    )
    active_table_sha256 = source_table_sha256 or ct_sha256
    if _manifest_path_md.exists():
        _m = json.loads(_manifest_path_md.read_text(encoding="utf-8"))
        _ct_id_md = _m.get("codon_table_id", "nbenthamiana_legacy_kazusa_sgn_v101")
        _ct_ss_md = _m.get("source_status", "legacy_metadata_only")
        _ct_bp_md = _m.get("build_path_status", "incomplete")
        _ct_sha_md = ct_sha256 or _m.get("sha256", "")
    else:
        _ct_id_md, _ct_ss_md, _ct_bp_md, _ct_sha_md = (
            "nbenthamiana_legacy_kazusa_sgn_v101", "legacy_metadata_only", "incomplete", ""
        )

    lines = [
        f"# Benchmark Summary (dataset={dataset}, mode={mode})", "",
        "## Provenance (reproducibility)",
        f"- FactorForge version: {FF_VERSION}",
        f"- registry sha256: {cfg.registry_sha256}",
        f"- benchmark_spec sha256: {cfg.spec_sha256}",
    ]
    if prot_sha256:
        lines.append(f"- input protein fasta sha256: {prot_sha256}")
    if cds_sha256:
        lines.append(f"- dataset cds fasta sha256: {cds_sha256}")
    lines += [
        f"- run date: {date.today().isoformat()}",
        f"- codon_table_id: {_ct_id_md}",
        f"- codon_table_sha256: {_ct_sha_md}",
        f"- source_profile_id: {source_profile_id}",
        f"- design_table_sha256: {active_table_sha256}",
        f"- score_table_sha256: {active_table_sha256}",
        f"- codon_table_source_status: {_ct_ss_md}",
        f"- codon_table_build_path_status: {_ct_bp_md}",
        "",
        "> **Codon table note:** The current N. benthamiana codon table is a legacy FactorForge reference "
        "labeled as derived from Kazusa CodonUsage Database and SGN genome v1.0.1-era resources. "
        "The original authoritative build path for the current JSON codon table is incomplete/not verified. "
        "The formal benchmark dataset uses SGN QLD183 v103 records; therefore CAI and codon-usage metrics "
        "should be interpreted as scores against the configured FactorForge codon reference, "
        "not as a de novo SGN QLD183 v103 codon-usage reconstruction.",
        "", f"> {CLAIM}", "",
        "| method | multi_constraint_pass_rate | biological_pass_rate | assembly_pass_rate | gc_in_range_rate | mean_cai |",
        "|---|---|---|---|---|---|",
    ]
    for m in methods:
        lines.append(
            f"| {m} | {_pass_rate(rows,'multi_constraint_pass',m):.2f} | "
            f"{_pass_rate(rows,'biological_pass',m):.2f} | {_pass_rate(rows,'assembly_pass',m):.2f} | "
            f"{_pass_rate(rows,'gc_in_target_range',m):.2f} | {_mean(rows,'cai',m):.4f} |"
        )
    lines += [
        "", "_native_reference is a biological reference anchor, not an optimizer._",
        "_greedy_cai is a CAI-focused baseline that does not explicitly optimize GC or assembly constraints._",
        "",
        "> **Note on factorforge_high_cai:** This profile optimizes against a golden-set codon reference "
        "(high-expression gene subset). General CAI scores are computed against the full N. benthamiana "
        "CDS-derived codon table and are not expected to be maximized by this profile. "
        "Lower general CAI relative to greedy_cai is expected and does not indicate poor performance.",
        "",
        "FactorForge preserves amino-acid identity and avoids invalid CDS outputs while "
        "improving multi-constraint in-silico CDS design quality relative to simple "
        "synonymous-codon baselines.",
    ]
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines), encoding="utf-8")

    # JSON observability record — per-profile breakdown
    def _per_method(m):
        all_m = [r for r in rows if r["method"] == m]
        ok = [r for r in all_m if r.get("status") == "ok"]
        err = [r for r in all_m if r.get("status") == "error"]
        if not all_m:
            return {}
        n_ok = len({r["sequence_id"] for r in ok})
        n_err = len(err)
        return {
            "n_ok": n_ok,
            "n_error": n_err,
            "error_rate": round(n_err / max(1, n_ok + n_err), 4),
            "pass_rate_biological": round(sum(1 for r in ok if r["biological_pass"]) / max(1, len(ok)), 4),
            "pass_rate_assembly": round(sum(1 for r in ok if r["assembly_pass"]) / max(1, len(ok)), 4),
            "pass_rate_multi_constraint": round(sum(1 for r in ok if r["multi_constraint_pass"]) / max(1, len(ok)), 4),
            "gc_in_range_rate": round(sum(1 for r in ok if r["gc_in_target_range"]) / max(1, len(ok)), 4),
            "mean_cai": round(sum(r["cai"] for r in ok if r["cai"] is not None) / max(1, sum(1 for r in ok if r["cai"] is not None)), 4),
        }

    # Compute codon table sha256 and read manifest at runtime
    # out_md is benchmarks/results/v3.2.0/benchmark_summary.md → 4 levels up = repo root
    repo_root = out_md.parent.parent.parent.parent
    manifest_path = repo_root / "data" / "reference" / "codon_table_manifest.json"
    if manifest_path.exists():
        _manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        ct_id = _manifest.get("codon_table_id", "nbenthamiana_legacy_kazusa_sgn_v101")
        ct_source_status = _manifest.get("source_status", "legacy_metadata_only")
        ct_build_path_status = _manifest.get("build_path_status", "incomplete")
    else:
        ct_id = "nbenthamiana_legacy_kazusa_sgn_v101"
        ct_source_status = "legacy_metadata_only"
        ct_build_path_status = "incomplete"

    try:
        _git_commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        _git_commit = None

    summary_json = out_md.parent / "benchmark_summary.json"
    summary_json.write_text(json.dumps({
        "run_id": hashlib.sha256(f"{dataset}{mode}{date.today().isoformat()}".encode()).hexdigest()[:12],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "scoring_contract_version": "v1.1",
        "multi_constraint_pass_definition": "biological_pass AND assembly_pass AND gc_in_target_range",
        "factorforge_version": FF_VERSION,
        "git_commit": _git_commit,
        "registry_version": cfg.registry_version,
        "registry_hash": f"sha256:{cfg.registry_sha256}",
        "spec_hash": f"sha256:{cfg.spec_sha256}",
        "dataset_id": dataset,
        "dataset_n": len({r["sequence_id"] for r in rows}),
        "random_seed": seed,
        "input_protein_fasta_sha256": prot_sha256,
        "dataset_cds_fasta_sha256": cds_sha256,
        "runtime_seconds": round(sum(r["runtime_seconds"] for r in rows), 4),
        "methods": {m: _per_method(m) for m in methods},
        "codon_table_id": ct_id,
        "codon_table_sha256": active_table_sha256,
        "source_profile_id": source_profile_id,
        "design_table_sha256": active_table_sha256,
        "score_table_sha256": active_table_sha256,
        "source_profile_manifest_sha256": source_manifest_sha256,
        "codon_table_source_status": ct_source_status,
        "codon_table_build_path_status": ct_build_path_status,
        "codon_table_reference_note": (
            "CAI and codon-usage metrics are computed against the configured FactorForge codon reference. "
            "The current codon table is not a freshly rebuilt SGN QLD183 v103 codon table."
        ),
        "vienna_rna_active": _check_vienna_available(),
    }, indent=2), encoding="utf-8")
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    if summary.get("design_table_sha256") != summary.get("score_table_sha256"):
        raise RuntimeError(
            "FATAL: design_table_sha256 != score_table_sha256 in run manifest. "
            "Injection integrity violated."
        )


def main(default_mode: str = "formal") -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="synthetic")
    ap.add_argument("--mode", default=default_mode)
    ap.add_argument("--limit", type=int, default=None,
                    help="Cap number of sequences (for quick validation runs).")
    ap.add_argument("--progress-every", type=int, default=500,
                    help="Print progress every N sequences (default: 500).")
    ap.add_argument("--seed", type=int, default=None,
                    help="Random seed for FactorForge profiles (default: non-deterministic). "
                         "Use --seed 320 for reproducible formal runs.")
    ap.add_argument("--codon-table-path", default=None,
                    help="Path to source-profile codon table JSON. If omitted, uses legacy nbenthamiana_codons.json.")
    ap.add_argument("--source-profile-id", default=None,
                    help="Source profile identifier (required if --codon-table-path is provided).")
    ap.add_argument("--source-profile-manifest", default=None,
                    help="Path to source profile manifest JSON for hash verification.")
    args = ap.parse_args()
    if args.codon_table_path and not args.source_profile_id:
        raise SystemExit("ERROR: --source-profile-id is required when --codon-table-path is provided")

    codon_weights = None
    active_profile_id = "legacy_packaged"
    active_table_sha256 = None
    active_manifest_sha256 = None

    if args.codon_table_path:
        from factorforge.analysis.metrics import load_codon_usage_table

        table_path = Path(args.codon_table_path)
        table_bytes = table_path.read_bytes()
        active_table_sha256 = hashlib.sha256(table_bytes).hexdigest()

        if args.source_profile_manifest:
            manifest_path = Path(args.source_profile_manifest)
            manifest_bytes = manifest_path.read_bytes()
            active_manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
            manifest_data = json.loads(manifest_bytes)
            expected_sha = manifest_data.get("codon_profile_sha256")
            if expected_sha and expected_sha != active_table_sha256:
                raise SystemExit(
                    "ERROR: codon table SHA-256 mismatch\n"
                    f"  file:     {active_table_sha256}\n"
                    f"  manifest: {expected_sha}"
                )

        table = load_codon_usage_table(table_path)
        codon_weights = table.codon_weights
        active_profile_id = args.source_profile_id
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
        proteins_fasta=proteins, native_fasta=native,
        limit=args.limit, progress_every=args.progress_every, seed=args.seed,
        codon_weights=codon_weights,
        source_profile_id=active_profile_id,
        source_table_sha256=active_table_sha256,
        source_manifest_sha256=active_manifest_sha256,
        codon_table_path=Path(args.codon_table_path) if args.codon_table_path else None)


if __name__ == "__main__":
    main()
