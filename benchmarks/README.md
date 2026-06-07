# FactorForge Benchmarks

In-silico benchmark foundation (v3.2.0). Evaluates CDS design quality only —
**not** protein expression/yield/folding (wet-lab validation required).

## Files
- `current_parameter_registry.yaml` — canonical parameter values + provenance
- `benchmark_spec.yaml` — run config (references the registry)
- `run_benchmark.py` — canonical runner
- `baselines/` — random_synonymous, most_frequent_codon, greedy_cai, native_reference
- `datasets/` — fetch script + provenance (raw data not committed)

## Quick smoke (no network)
`python -m pytest tests/test_benchmark_smoke.py -v`

## Formal run (fetches real data)
`python benchmarks/datasets/fetch_dataset.py --url <verified-source>`
`python benchmarks/run_benchmark.py --dataset nbenthamiana_full --mode formal`
Results: `benchmarks/results/v3.2.0/`

## Result freezing
Results under `results/v3.2.0/` are frozen after release. If benchmark logic or
dataset changes, create a new versioned results directory (e.g. `results/v3.2.1/`).
native_reference is a biological reference anchor, not an optimizer.
