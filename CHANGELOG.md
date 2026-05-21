# Changelog

All notable changes to FactorForge are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased] — v3-alpha

### Added
- **Metrics Engine** — CAI, GC (global/local/first-region), amino acid identity, internal stop count, homopolymer/repeat/forbidden motif/invalid codon detection (`src/factorforge/ml/metrics.py`)
- **Sequence Validator** — structured hard-fail contract (AA identity, stops, invalid codons, GC, length) returning machine-readable dict (`src/factorforge/utils/validation.py`)
- **Pareto Feasibility Analyzer** — computes achievable CAI/GC ranges under synonymous codon constraints; distinguishes infeasible objectives from model failure (`src/factorforge/ml/feasibility.py`)
- **v2 Adapter** — formal baseline wrapper, returns structured dict with metrics and validator result (`src/factorforge/engines/v3/inference/v2_adapter.py`)
- **Synonym Mask** — per-position boolean mask restricting decoding to valid synonymous codons (`src/factorforge/engines/v3/synonym_mask.py`)
- **Constrained Decoder** — synonym-constrained autoregressive decoding + post-validation + v2 fallback on failure (`src/factorforge/engines/v3/inference/constrained_decoder.py`)
- **v2 Pseudo-label Generator** — generates Run 4 training targets from v2 optimizer (`scripts/1_data_preparation/generate_v2_pseudolabels.py`)
- **Run 4 Config** — `configs/v3_training_config_alpha_run1.yml`: CE-to-pseudo-label + bounded GC penalty + expected log CAI reward + synonym mask required
- **Benchmark Panel** — 5 synthetic proteins in `tests/fixtures/benchmark_proteins.py`
- **Diagnostic & benchmark scripts** — `diagnose_v2_v3_metrics.py`, `benchmark_candidates.py`, `check_constraint_feasibility.py`, `run_v2_baseline_benchmark.py`
- **Docs** — `docs/reports/v3_audit_report.md`, `v3_loss_debug_report.md`, `docs/design/v3_alpha_product_boundary.md`, `docs/model_cards/factorforge_v3_alpha_model_card.md`

### Changed
- GC loss upgraded: single-target absolute penalty → bounded GC range penalty (`gc_low=0.40, gc_high=0.55`)
- Training loss extended with differentiable expected log CAI reward
- Synonym mask applied in expected GC/log CAI loss computation
- `train_v3_esm2_bart.py` — bounded GC + expected log CAI + synonym mask support
- `dataset.py` — amino-acid-aware synonym masks aligned to decoder labels per batch

### Key Findings
- GFP CAI ≥ 0.92 + GC 41–44% simultaneously is **infeasible** under current N. benthamiana codon weights (Pareto max CAI at GC 41–44%: 0.884; at GC 40–55%: 0.969)
- GC collapse (32–33%) root cause: native CDS CE objective mismatch + single-target GC penalty + no expected CAI reward + no synonym mask + special-token probability leakage
- Synonym mask + constrained decoding now prevent non-synonymous codon output and internal stop codons

### Run 2 / Run 3 Results (N. benthamiana, GFP)
| Metric | v2 Baseline | v3 Run 1 | v3 Run 2 | v3 Run 3 |
|--------|------------|---------|---------|---------|
| CAI | 0.909 | 0.974 | 0.717 | 0.718 |
| GC% | 56.4% | 60.9% | 32.6% | 30.9% |
| Eval loss | — | 6.66 | 2.505 | 2.734 |
| Training data | — | 917 | 38,754 | 38,754 |

### Training Results

**alpha_run1** (v3-alpha, *N. benthamiana*, pseudo-label high_cai, 34,878 sequences):
| Metric | v2 baseline | v3 alpha_run1 |
|--------|------------|--------------|
| CAI mean | 0.733 | 0.733 |
| GC% mean | 33.23% | 33.24% |
| Validator pass rate | — | 100% |
| AA identity | — | 100% |

Outcome: **GC target miss** — decoded GC remained below 40% target (mean 33.24%). Root cause: `gc_weight=0.2` too weak relative to CE loss. See `experiments/results/alpha_run1/`.

**alpha_run2** (v3-alpha, *N. benthamiana*, pseudo-label mixed, ~30,000 sequences):
| Metric | v2 baseline | v3 alpha_run2 |
|--------|------------|--------------|
| CAI mean | 0.8025 | 0.8025 |
| GC% mean | 42.54% | 42.54% |
| GC% range | 40.36–53.81% | 41.99–53.81% |
| Validator pass rate | — | 100% |
| AA identity | — | 100% |
| Eval sequences | — | 3,876 |

Outcome: **✅ GC recovery** — mixed teacher distribution + `gc_weight=2.0` resolved GC collapse. See `experiments/results/alpha_run2/`.

### Tests
- 299 → 333 passing (34 new tests: metrics, validation, feasibility, synonym mask, constrained decoder, loss gradient flow)

---

## [1.0.0] — 2026-05-18

Initial public release of FactorForge.

### Added
- **ESM2 + BART decoder architecture** — per-token protein embeddings (esm2_t6_8M_UR50D, 320-dim) → autoregressive codon generation
- **v2 rule-based optimizer** — CAI-maximizing baseline (RuleBasedOptimizer v2.5.2)
- **SGN data pipeline** — fetch, CAI-filter, and split N. benthamiana CDS from Sol Genomics Network v2.6.1
- **GC penalty loss** — differentiable GC content penalty targeting N. benthamiana spec (41–44%)
- **Biological regression tests** — GFP CAI ≥ 0.80, GC 41–44%, no PolyA signals
- **Codon metrics module** — CAI, GC%, GC3%, ENC, homopolymer detection, rare codon runs
- **Baseline comparison script** — v2 vs v3 side-by-side (CAI, GC%, ENC)

### Architecture
- Encoder: ESM2 esm2_t6_8M_UR50D (frozen, 320-dim per-token)
- Decoder: BART (6 layers, d_model=256, vocab_size=69)
- Parameters: 4,206,336
- Host: Nicotiana benthamiana (primary)

### Run 1 Results (N. benthamiana, GFP)
| Metric | v2 Baseline | v3 Run 1 |
|--------|------------|---------|
| CAI | 0.909 | 0.974 |
| GC% | 56.4% | 60.9% |
| Training data | — | 917 sequences |
| Eval loss | — | 6.66 (overfitting) |

v3-alpha training in progress. alpha_run1 (GC target miss) and alpha_run2 (✅ GC recovery) completed — see [Unreleased].

---

[Unreleased]: https://github.com/eijex/factorforge/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/eijex/factorforge/releases/tag/v1.0.0
