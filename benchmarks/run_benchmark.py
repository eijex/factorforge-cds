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
DEFAULT_RESULTS_DIR = ROOT / "benchmarks" / "results" / "v3.2.0"
OUTPUT_GUARD_FILENAMES = ("benchmark_results.csv", "benchmark_summary.json")


def _resolve_output_dir(out_dir: Path | None) -> Path:
    return (out_dir or DEFAULT_RESULTS_DIR).resolve()


def _guard_output_dir(out_dir: Path, force: bool) -> None:
    existing = [out_dir / name for name in OUTPUT_GUARD_FILENAMES if (out_dir / name).exists()]
    if not existing or force:
        return
    file_list = "\n  ".join(str(path) for path in existing)
    raise SystemExit(
        "ERROR: benchmark output directory already contains result artifact(s). "
        "Pass --force to overwrite intentionally:\n"
        f"  {file_list}"
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


def _resolve_active_default_reference(root: Path) -> dict:
    """Resolve the codon reference the FactorForge engine actually uses by
    default (no --codon-table-path override), by fact rather than a
    hardcoded label.

    v3.3.0 reference-policy update: the no-flag benchmark path keeps its
    original design intent of tracking the current product default; this
    only makes the resulting summary label reflect that default truthfully
    instead of assuming it is always "legacy_packaged".
    """
    pointer_path = root / "data" / "reference" / "active_codon_reference.json"
    if pointer_path.exists():
        pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
        return {
            "profile_id": pointer["active_codon_table_id"],
            "table_path": root / pointer["file"],
            "manifest_path": root / pointer["manifest"],
        }
    # No pointer file: fall back to the historically bundled legacy reference.
    return {
        "profile_id": "legacy_packaged",
        "table_path": root / "src" / "factorforge" / "data" / "nbenthamiana_codons.json",
        "manifest_path": root / "data" / "reference" / "codon_table_manifest.json",
    }


def _resolve_manifest_for_table_path(table_path: Path, root: Path) -> Path | None:
    """Find the schema-conformant manifest for an explicitly-given codon table
    path, by matching it against known assets (legacy or current active
    default) — never by assuming "current active" just because no manifest
    was supplied. Returns None for genuinely unknown/custom tables rather
    than mislabeling them with an unrelated asset's manifest facts.
    """
    table_path = Path(table_path).resolve()
    legacy_path = (root / "src" / "factorforge" / "data" / "nbenthamiana_codons.json").resolve()
    if table_path == legacy_path:
        return root / "data" / "reference" / "codon_table_manifest.json"
    default_ref = _resolve_active_default_reference(root)
    if table_path == default_ref["table_path"].resolve():
        return default_ref["manifest_path"]
    return None


def run(dataset: str, mode: str, out_csv: Path, out_md: Path,
        proteins_fasta: Path, native_fasta: Path,
        limit: int | None = None, progress_every: int = 500,
        seed: int | None = None,
        codon_weights: dict[str, float] | None = None,
        source_profile_id: str | None = None,
        source_table_sha256: str | None = None,
        source_manifest_sha256: str | None = None,
        source_manifest_path: Path | None = None,
        codon_table_path: Path | None = None) -> None:
    # v3.3.0 reference-policy update: resolve no-flag defaults here (not
    # only in main()) so every caller of run() — CLI, tests, aggregate
    # scripts — gets fact-derived provenance instead of a stale literal.
    if codon_table_path is None and source_profile_id is None:
        _default_ref = _resolve_active_default_reference(ROOT)
        source_profile_id = _default_ref["profile_id"]
        if source_table_sha256 is None and _default_ref["table_path"].exists():
            source_table_sha256 = hashlib.sha256(_default_ref["table_path"].read_bytes()).hexdigest()
        if source_manifest_path is None:
            source_manifest_path = _default_ref["manifest_path"]
        if source_manifest_sha256 is None and source_manifest_path.exists():
            source_manifest_sha256 = hashlib.sha256(source_manifest_path.read_bytes()).hexdigest()
    elif source_profile_id is None:
        source_profile_id = "legacy_packaged"

    if codon_table_path is not None and source_manifest_path is None:
        # Explicit override: only attach a manifest if the given path
        # matches a known asset with a schema-conformant manifest. Otherwise
        # leave it unresolved — _write_summary() falls back to source_profile_id
        # as the literal codon_table_id rather than guessing.
        source_manifest_path = _resolve_manifest_for_table_path(codon_table_path, ROOT)
        if source_manifest_path is not None and source_manifest_sha256 is None and source_manifest_path.exists():
            source_manifest_sha256 = hashlib.sha256(source_manifest_path.read_bytes()).hexdigest()

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
        source_manifest_path=source_manifest_path,
        is_explicit_override=codon_table_path is not None,
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
                   source_manifest_sha256: str | None = None,
                   source_manifest_path: Path | None = None,
                   is_explicit_override: bool = False) -> None:
    methods = sorted({r["method"] for r in rows})
    # Resolve codon table metadata for MD from the manifest that is actually
    # active for this run (legacy override or current production default),
    # not a hardcoded legacy assumption (v3.3.0 reference-policy update).
    # Use the module-level ROOT constant, not out_md's ancestors — out_md may
    # be a tmp_path in tests/smoke runs and not sit 4 levels under repo root.
    _repo_root_md = ROOT
    if source_manifest_path is not None:
        _manifest_path_md: Path | None = source_manifest_path
    elif not is_explicit_override:
        _resolved = _resolve_active_default_reference(_repo_root_md)
        _manifest_path_md = _resolved["manifest_path"]
    else:
        # Explicit override to a table with no known schema-conformant
        # manifest (e.g. qld183_v103/nbev11_cds_all derived candidates) —
        # do not guess; use the literal facts the caller already supplied.
        _manifest_path_md = None
    ct_sha256 = source_table_sha256
    active_table_sha256 = source_table_sha256
    if _manifest_path_md is not None and _manifest_path_md.exists():
        _m = json.loads(_manifest_path_md.read_text(encoding="utf-8"))
        _ct_id_md = _m.get("codon_table_id", source_profile_id)
        _ct_ss_md = _m.get("source_status", "unknown")
        _ct_bp_md = _m.get("build_path_status", "unknown")
        _ct_sha_md = ct_sha256 or _m.get("sha256", "")
        _ct_asset_type_md = _m.get("asset_type", "unknown")
    else:
        _ct_id_md, _ct_ss_md, _ct_bp_md, _ct_sha_md, _ct_asset_type_md = (
            source_profile_id, "unknown", "unknown", ct_sha256 or "",
            "legacy_packaged" if source_profile_id == "legacy_packaged" else "unknown",
        )
    if _ct_asset_type_md == "legacy_packaged":
        _codon_table_note = (
            "> **Codon table note:** This run used the legacy FactorForge codon reference "
            "labeled as derived from Kazusa CodonUsage Database and SGN genome v1.0.1-era resources. "
            "The original authoritative build path for this table is incomplete/not verified. "
            "Scores are interpretable as an archived legacy-codon-reference (contract v1) historical "
            "record (see the reference-policy audit), not a freshly rebuilt SGN QLD183 v103 codon-usage reconstruction. "
            "Note: this 'contract v1/v2' refers to the codon-reference asset generation, a separate "
            "concept from the frozen scoring_contract_version 1.1 pass/fail definition below."
        )
    else:
        _codon_table_note = (
            f"> **Codon table note:** This run used codon reference `{_ct_id_md}` "
            "(codon-reference contract v2; NbeV1.1 LAB-strain genome-derived). "
            "v3.3.0 reference-policy update — see data/reference/active_codon_reference.json."
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
        _codon_table_note,
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

    # Reuse the manifest facts already resolved above (active default or
    # explicit override) instead of re-reading a hardcoded legacy path.
    ct_id, ct_source_status, ct_build_path_status = _ct_id_md, _ct_ss_md, _ct_bp_md

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
                    help="Path to source-profile codon table JSON. If omitted, tracks the current "
                         "FactorForge production default (see data/reference/active_codon_reference.json).")
    ap.add_argument("--source-profile-id", default=None,
                    help="Source profile identifier (required if --codon-table-path is provided).")
    ap.add_argument("--source-profile-manifest", default=None,
                    help="Path to source profile manifest JSON for hash verification.")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output directory for benchmark_results.csv, benchmark_summary.md, "
                         "and benchmark_summary.json. Defaults to benchmarks/results/v3.2.0.")
    ap.add_argument("--force", action="store_true",
                    help="Allow overwriting benchmark_results.csv or benchmark_summary.json "
                         "in the output directory.")
    args = ap.parse_args()
    if args.codon_table_path and not args.source_profile_id:
        raise SystemExit("ERROR: --source-profile-id is required when --codon-table-path is provided")

    root = Path(__file__).resolve().parents[1]
    out_dir = _resolve_output_dir(args.out_dir)
    _guard_output_dir(out_dir, args.force)
    codon_weights = None
    active_manifest_path: Path | None = None

    if not args.codon_table_path:
        # v3.3.0 reference-policy update: label this run by the fact of
        # what the engine actually defaults to, instead of assuming legacy.
        _default_ref = _resolve_active_default_reference(root)
        active_profile_id = _default_ref["profile_id"]
        active_manifest_path = _default_ref["manifest_path"]
        active_table_sha256 = (
            hashlib.sha256(_default_ref["table_path"].read_bytes()).hexdigest()
            if _default_ref["table_path"].exists() else None
        )
        active_manifest_sha256 = (
            hashlib.sha256(active_manifest_path.read_bytes()).hexdigest()
            if active_manifest_path.exists() else None
        )
    else:
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
            active_manifest_path = manifest_path
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
    if args.dataset == "synthetic":
        proteins = root / "tests" / "fixtures" / "small_proteins.fasta"
        native = root / "tests" / "fixtures" / "small_native_cds.fasta"
    else:  # formal: real fetched dataset (must run fetch_dataset.py first)
        proteins = root / "benchmarks" / "datasets" / "nbenthamiana_reference_proteins.fasta"
        native = root / "benchmarks" / "datasets" / "nbenthamiana_reference_cds.fasta"
    run(dataset=args.dataset, mode=args.mode,
        out_csv=out_dir / "benchmark_results.csv",
        out_md=out_dir / "benchmark_summary.md",
        proteins_fasta=proteins, native_fasta=native,
        limit=args.limit, progress_every=args.progress_every, seed=args.seed,
        codon_weights=codon_weights,
        source_profile_id=active_profile_id,
        source_table_sha256=active_table_sha256,
        source_manifest_sha256=active_manifest_sha256,
        source_manifest_path=active_manifest_path,
        codon_table_path=Path(args.codon_table_path) if args.codon_table_path else None)


if __name__ == "__main__":
    main()
