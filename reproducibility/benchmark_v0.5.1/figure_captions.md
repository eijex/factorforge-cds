# FactorForge Benchmark Figure Captions

> Source: `reproducibility/benchmark_v0.5.1/`
> Version: FactorForge v3.2.0 | N=49,257 | seed=320 | scoring_contract v1.1
> All metrics are computational sequence-design constraint satisfaction only.
> No wet-lab expression or yield data is included or implied.

---

## Figure 2 — Multi-Constraint Pass Rate by Method

**File**: `figures/figure2_multiconstraint_pass_rate.{png,svg}`

**Caption**:
Multi-constraint pass rates (biological_pass AND assembly_pass AND gc_in_target_range) across optimization methods. FactorForge v3.2.0, N. benthamiana SGN CDS, N=49,257, balanced/assembly_friendly profile, seed=320, scoring_contract v1.1. Random synonymous: mean of 3 replicates (seeds 320/321/322). Metrics represent computational sequence-design constraint satisfaction only and do not measure wet-lab expression level or yield.

---

## Figure 3 — Method Comparison Heatmap

**File**: `figures/figure3_benchmark_tradeoff_heatmap.{png,svg}`

**Caption**:
Multi-constraint metric heatmap across optimization methods (FactorForge v3.2.0, N. benthamiana SGN CDS, N=49,257, seed=320, scoring_contract v1.1). Rows: methods. Columns: biological_pass, assembly_pass, gc_in_target_range, multi_constraint_pass. All metrics are in-silico; no wet-lab validation has been performed.

---

## Figure 4 — Failure Attribution

**File**: `figures/figure_failure_attribution.png`

**Caption**:
Constraint-level failure attribution for FactorForge v3.2.0 (N. benthamiana SGN CDS, N=49,257, seed=320, scoring_contract v1.1). Categories: gc_out_of_range, assembly_fail_only, both_fail, biological_fail, pass. Counts are computational only and do not reflect wet-lab outcomes.

---

## Figure 5 — CAI Distribution

**File**: `figures/figure_cai_distribution.png`

**Caption**:
Distribution of Codon Adaptation Index (CAI) across optimization methods (FactorForge v3.2.0, N. benthamiana SGN CDS, N=49,257, seed=320). Reference codon usage: N. benthamiana (Kazusa + SGN v1.0.1, see `data/reference/codon_table_manifest.json`). CAI is a sequence-level metric; it does not predict in-vivo expression level.

---

## Figure 6 — GC Content Distribution

**File**: `figures/figure_gc_distribution.png`

**Caption**:
Distribution of CDS-level GC content across optimization methods (FactorForge v3.2.0, N. benthamiana SGN CDS, N=49,257, seed=320). Target range: 55–65% (registry: `gc_range_nbenthamiana_global`). GC compliance is a computational constraint; it does not guarantee expression or cloning success.

---

## Standard Footer (all figures)

> Metrics represent computational constraint satisfaction only and do not predict wet-lab expression, yield, or cloning success.
> Source data: `benchmark_summary.frozen.json` (FactorForge v3.2.0, scoring_contract v1.1).
> Reproducibility: see `MANIFEST.json`.
