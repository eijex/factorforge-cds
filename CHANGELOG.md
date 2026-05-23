# Changelog

All notable changes to FactorForge are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added

- Custom restriction site detection and synonymous CDS cleanup for API and web UI requests.
- Rare codon run scanner in the v2 rule engine, included in full scan mode and excluded from fast scan mode.
- Greedy CpG/TpA dinucleotide reduction pass in the v2 rule engine and optimization pipeline.

---

## [3.0.0] — 2026-05-23

First official release of FactorForge.

### Distribution
- **PyPI** — `pip install factorforge-cds` ([pypi.org/project/factorforge-cds](https://pypi.org/project/factorforge-cds/))
- **Docker** — `docker pull ghcr.io/eijex/factorforge-cds:latest` for local deployment

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

### Tests

- 333 passing (metrics, validation, feasibility, synonym mask, constrained decoder, loss gradient flow)

---

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.0.0...HEAD
[3.0.0]: https://github.com/eijex/factorforge-cds/releases/tag/v3.0.0
