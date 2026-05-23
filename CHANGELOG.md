# Changelog

All notable changes to FactorForge are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [3.0.0] — 2026-05-23

First official release of FactorForge.

### Added

- **v2 Rule-based Optimizer** — CAI-maximizing baseline with profiles: `balanced`, `high_cai`, `gc_target`, `viral_delivery` (`src/factorforge/engines/v2/`)
- **DP Feasibility Engine** — constraint-based dynamic programming; computes achievable CAI/GC ranges under synonymous codon constraints (`src/factorforge/ml/feasibility.py`)
- **Metrics Engine** — CAI, GC (global/local/first-region), amino acid identity, internal stop count, homopolymer/repeat/forbidden motif/invalid codon detection (`src/factorforge/ml/metrics.py`)
- **Sequence Validator** — structured hard-fail contract (AA identity, stops, invalid codons, GC, length) returning machine-readable dict (`src/factorforge/utils/validation.py`)
- **Synonym Mask** — per-position boolean mask restricting decoding to valid synonymous codons (`src/factorforge/engines/v3/synonym_mask.py`)
- **Constrained Decoder** — synonym-constrained autoregressive decoding + post-validation + v2 fallback on failure (`src/factorforge/engines/v3/inference/constrained_decoder.py`)
- **v2 Adapter** — formal baseline wrapper, returns structured dict with metrics and validator result (`src/factorforge/engines/v3/inference/v2_adapter.py`)
- **SGN data pipeline** — fetch, CAI-filter, and split *N. benthamiana* CDS from Sol Genomics Network v2.6.1
- **Benchmark Panel** — 5 synthetic proteins in `tests/fixtures/benchmark_proteins.py`
- **Diagnostic & benchmark scripts** — `diagnose_v2_v3_metrics.py`, `benchmark_candidates.py`, `check_constraint_feasibility.py`, `run_v2_baseline_benchmark.py`
- **CLI** — `factorforge optimize` / `factorforge list-engines` with FASTA and GenBank output
- **Docs** — `docs/reports/v3_audit_report.md`, `docs/design/v3_alpha_product_boundary.md`, `docs/model_cards/factorforge_v3_alpha_model_card.md`

### Architecture

- Encoder: ESM2 esm2_t6_8M_UR50D (frozen, 320-dim per-token)
- Decoder: BART (6 layers, d_model=256, vocab_size=69)
- Parameters: 4,206,336
- Host: *Nicotiana benthamiana* (primary)

### Key Findings

- GFP CAI ≥ 0.92 + GC 41–44% simultaneously is **infeasible** under current *N. benthamiana* codon weights (Pareto max CAI at GC 41–44%: 0.884; at GC 40–55%: 0.969)
- GC collapse (32–33%) root cause: native CDS CE objective mismatch + single-target GC penalty + no expected CAI reward + no synonym mask + special-token probability leakage
- Synonym mask + constrained decoding now prevent non-synonymous codon output and internal stop codons

### ML Training Results (v3-alpha, *N. benthamiana*)

**Run 1–3** (GFP, early experiments):

| Metric | v2 Baseline | v3 Run 1 | v3 Run 2 | v3 Run 3 |
|--------|-------------|---------|---------|---------|
| CAI | 0.909 | 0.974 | 0.717 | 0.718 |
| GC% | 56.4% | 60.9% | 32.6% | 30.9% |
| Eval loss | — | 6.66 | 2.505 | 2.734 |
| Training data | — | 917 | 38,754 | 38,754 |

**alpha_run1** (pseudo-label high_cai, 34,878 sequences):

| Metric | v2 baseline | v3 alpha_run1 |
|--------|-------------|--------------|
| CAI mean | 0.733 | 0.733 |
| GC% mean | 33.23% | 33.24% |
| Validator pass rate | — | 100% |
| AA identity | — | 100% |

Outcome: **GC target miss** — decoded GC remained below 40% target. Root cause: `gc_weight=0.2` too weak relative to CE loss.

**alpha_run2** (pseudo-label mixed, ~30,000 sequences):

| Metric | v2 baseline | v3 alpha_run2 |
|--------|-------------|--------------|
| CAI mean | 0.8025 | 0.8025 |
| GC% mean | 42.54% | 42.54% |
| GC% range | 40.36–53.81% | 41.99–53.81% |
| Validator pass rate | — | 100% |
| AA identity | — | 100% |
| Eval sequences | — | 3,876 |

Outcome: **✅ GC recovery** — mixed teacher distribution + `gc_weight=2.0` resolved GC collapse.

### Tests

- 333 passing (metrics, validation, feasibility, synonym mask, constrained decoder, loss gradient flow)

---

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.0.0...HEAD
[3.0.0]: https://github.com/eijex/factorforge-cds/releases/tag/v3.0.0
