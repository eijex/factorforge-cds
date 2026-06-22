# FactorForge Benchmark Evidence Pack v0.5.1

> **Correction Notice (2026-06-13)**: `multi_constraint_pass` definition corrected in scoring_contract v1.1.
> Previous definition (v1.0): `biological_pass AND assembly_pass` — GC target compliance was omitted.
> Corrected definition: `biological_pass AND assembly_pass AND gc_in_target_range`.
> All artifacts in this directory reflect the corrected definition. L3/L4 ablation values were previously inflated.
> Zenodo `benchmark_results.csv` v1 (DOI: 10.5281/zenodo.20640931) is superseded by v2 (DOI: 10.5281/zenodo.20676276, scoring_contract v1.1 corrected).

## Artifact Lineage

| Field | Value |
|-------|-------|
| FactorForge version | 3.2.0 |
| Registry version | 3.2.0 |
| Registry hash | sha256:8838408157e3e2913c2eb2de2168b369d89cd8ddff9a8c85c4bbced29532cc06 |
| Codon table ID | nbenthamiana_legacy_kazusa_sgn_v101 |
| Codon table SHA256 | ddbd0a41da88109a709bca0304581e29bd0a756e4db1c51809d5002e9b2d5e8c |
| Dataset N | 49257 |
| Random seed | 320 |
| Evidence-pack manifest git commit | 6409108921d538396b67c7edd3434b43e547f666 |
| Corrected formal summary git commit | e57341c |
| Corrected formal summary run ID | 3451bd347dcd |
| benchmark_summary.json SHA256 | f8f5aa8a2b61b04139bb1d77547e6e3887f46ceca350cd834dc306d9a111956a |
| scoring_contract_version | v1.1 (multi_constraint_pass = bio AND assembly AND gc_in_target_range) |
| Methods | random_synonymous, greedy_cai, native_reference, factorforge_balanced, factorforge_gc_target, factorforge_high_cai, factorforge_assembly_friendly |

The ablation summary retains its historical source identifiers
(`source_formal_run_id=c8641f5bbc32`,
`source_formal_summary_sha256=51aecac78ca3945aea977a62ac436b06400ddacc94e12bbe8486111c56728c21`)
for lineage. The frozen formal summary in this package is the corrected
scoring-contract v1.1 aggregate identified above; the two hashes should not be
treated as interchangeable.

## Evidence Boundary

All metrics in this package are **computational design constraint outcomes** derived from
FactorForge's deterministic constraint-aware CDS design harness. They measure:

- Amino-acid identity preservation
- Absence of internal stop codons and invalid codons
- GC content within target range (configured for *N. benthamiana*)
- Absence of forbidden Type IIS restriction enzyme recognition sites

**These metrics do not constitute wet-lab expression validation, yield prediction,
synthesis acceptance likelihood, or cloning-success claims.**

## Contents

| File/Folder | Description |
|-------------|-------------|
| `figures/figure2_multiconstraint_pass_rate.{png,svg}` | Evidence Pack Figure 2 — multi-constraint pass rate by method |
| `figures/figure3_benchmark_tradeoff_heatmap.{png,svg}` | Evidence Pack Figure 3 — metric heatmap |
| `tables/table3_benchmark_summary.md` | Evidence Pack Table 3 |
| `data/benchmark_summary.frozen.json` | Checksum-pinned copy of source benchmark_summary.json |
| `data/benchmark_summary_sha256.json` | SHA256 of source summary |
| `data/input_audit.json` | CSV row-count and method-list audit |
| `data/failure_attribution_summary.json` | Failure attribution per method |
| `data/constraint_overlap_matrix.csv` | GC × Type IIS co-failure matrix |
| `data/bootstrap_ci.json` | Sequence-level paired bootstrap CI |

## Legacy

The earlier v0.4 reproducibility package (parent repo) contains the v0.4 smoke-test artifacts.
Do not delete. The v0.5.1 package supersedes it for evidence-pack purposes but does not replace it.

## P2 Generation-Layer Ablation

### What the ablation measures

The ablation isolates the contribution of each constraint layer to the overall
multi-constraint pass rate, using the same `score_cds()` function applied uniformly
across all layers. Layers are defined by which constraints are active:

> **Note**: These values reflect scoring_contract v1.1 (multi_constraint_pass includes gc_in_target_range).
> Previous v1.0 values (omitting GC) were L3=89.0%, L4=88.6% — mathematically impossible given GC in range rates.

| Layer | Constraints active | Method name | Multi-constraint pass rate | GC in range rate |
|-------|--------------------|-------------|---------------------------|-----------------|
| L0 | None (random) | `random_synonymous` | 0.5% | 1.3% |
| L1 | CAI only (greedy) | `greedy_cai` | 21.8% | 72.9% |
| L2 | CAI + GC target | `ablation_cai_gc` | 26.0% | 99.8% |
| L3 | CAI + TypeIIS avoidance | `ablation_cai_type_iis` | 3.5% | 3.7% |
| L4 | CAI + GC + TypeIIS | `ablation_cai_gc_type_iis` | 5.6% | 5.8% |
| L5 | FactorForge full (assembly-friendly) | `factorforge_assembly_friendly` | 63.4% | 97.1% |

Dataset: N=49,257 *N. benthamiana* reference proteins. Seed: 320.

### Key finding (scoring_contract v1.1)

**GC target compliance is the dominant constraint**: multi_constraint_pass requires all three
conditions (biological_pass, assembly_pass, AND gc_in_target_range). The corrected staircase shows:

- L0 (random): 0.5% — baseline, virtually no sequences pass all three constraints by chance
- L1 (CAI only): 21.8% — CAI optimization increases GC in range rate (72.9%), assembly pass rate improves too
- L2 (CAI+GC): 26.0% — adding explicit GC steering achieves 99.8% GC compliance, but TypeIIS sites remain (~72% assembly pass) limiting multi-constraint pass
- L3 (CAI+TypeIIS): 3.5% — TypeIIS avoidance without GC steering: assembly pass is high (~89%) but GC in range is only 3.7%, so multi-constraint pass collapses
- L4 (CAI+GC+TypeIIS): 5.6% — both GC and TypeIIS attempted simultaneously; marginal improvement over L3 with GC targeting, but max_attempts=50 insufficient to satisfy all constraints simultaneously in most sequences
- L5 (FactorForge full): 63.4% — the production optimizer achieves the highest multi-constraint pass rate by combining balanced base composition (implicit GC steering) with TypeIIS avoidance

The key insight: **TypeIIS avoidance AND GC compliance are jointly difficult to satisfy**. The production optimizer (L5) achieves ~63% by using a different algorithmic approach (balanced base composition) rather than direct GC minimization + TypeIIS avoidance.

### Evidence boundary

These results are **in-silico constraint ablation outcomes** on a single reference proteome
(*N. benthamiana*). They quantify which constraint layers drive the computational
multi-constraint pass rate metric. They do not constitute:

- wet-lab expression or assembly validation
- synthesis success rate prediction
- claims about other host organisms

### Artifacts

| File | Description |
|------|-------------|
| `benchmarks/results/v3.2.0/ablation/ablation_summary.json` | Layer-by-layer metric summary (git-tracked) |
| `benchmarks/results/v3.2.0/ablation/figures/figure_ablation_pass_rate.{png,svg}` | Pass rate by layer bar chart |
| `benchmarks/results/v3.2.0/ablation/figures/figure_ablation_tradeoff_heatmap.{png,svg}` | Metric heatmap across layers |
| `benchmarks/results/v3.2.0/ablation/ablation_results.csv` | Full row-level results (74 MB — gitignored; regenerate via `run_ablation.py`) |

### Reproducing the ablation

```bash
# From repo root:
python -m benchmarks.ablation.run_ablation
```

`ablation_results.csv` is excluded from git (74 MB). The summary JSON and figures
are the canonical reproducibility artifacts for this layer.

---

## Reproducing the Formal Benchmark

The raw SGN archive is not committed. Place the downloaded
`NbQld183.v103.gff3.CDS.fasta.gz` file under `benchmarks/datasets/`, then run:

```bash
# From repo root:
python benchmarks/datasets/fetch_dataset.py \
  --file benchmarks/datasets/NbQld183.v103.gff3.CDS.fasta.gz
python -m benchmarks.run_benchmark \
  --dataset nbenthamiana_full --mode formal --seed 320
```

The formal benchmark scores amino-acid identity/internal-stop validity, global
GC target compliance, and forbidden Type IIS site absence. PolyA motifs, local
GC windows, rare-codon patterns, repeats, and homopolymers are reported by other
FactorForge surfaces but are not part of this formal benchmark scoring contract.

---

## Generating Figures

```bash
# From repo root:
python reproducibility/benchmark_v0.5.1/scripts/figures/make_benchmark_figures.py
```

Requires only `benchmarks/results/v3.2.0/benchmark_summary.json`. The 47 MB
`benchmark_results.csv` is not required for Figures 2–3.
If `benchmark_results.csv` is unavailable locally, reference it via Zenodo DOI
`10.5281/zenodo.20676276` (scoring_contract v1.1, corrected).
