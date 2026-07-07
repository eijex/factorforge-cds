# Benchmark Workflow

FactorForge benchmark runs are in-silico evidence artifacts. They do not
demonstrate expression, yield, folding, or wet-lab performance.

## 1. Smoke

Use synthetic or tiny fixture inputs to verify runner wiring, output schema,
hash fields, and local environment health. Smoke runs must write to a
throwaway local directory, not to `benchmarks/results/`, and must never create
ValidationHub records.

Example:

```bash
python benchmarks/run_benchmark.py --dataset synthetic --mode smoke \
  --limit 2 --out-dir tmp/benchmark-smoke --force
```

## 2. Preflight

Before a formal run, verify that the dataset files, codon-reference assets,
manifest files, intended `--out-dir`, git commit, and Python environment are
the intended source of truth. Preflight produces a local note or checklist only.
It does not write ValidationHub records.

## 3. Formal

A formal run uses the real benchmark dataset, an explicit seed, and a new
versioned output directory. The runner refuses to overwrite
`benchmark_results.csv` or `benchmark_summary.json` unless `--force` is passed.
Use `--force` only for an intentional rerun into the same directory.

Example:

```bash
python benchmarks/run_benchmark.py --dataset nbenthamiana_full --mode formal \
  --seed 320 --out-dir benchmarks/results/vNEXT
```

The formal output directory is the artifact source for CSV, Markdown summary,
JSON summary, manifests, and digest checks.

## 4. Evidence Register

Evidence registration happens only after a completed formal run. A human must
decide whether to create a ValidationHub `benchmark_evidence_reference` record.
This step is always human-gated and is never automatic.

The current minimal ValidationHub reference schema is unchanged. It records
only:

- `benchmark_run_id`
- `artifact_uri`
- `manifest_digest`

Do not create richer benchmark evidence fields here. Broader claim-boundary
metadata belongs in a future, deliberately scoped claim-evidence package.
