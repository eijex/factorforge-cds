# Changelog

All notable changes to FactorForge are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Versioning Policy

| Bump | When to use |
|------|-------------|
| **Major** (`X.0.0`) | Breaking API change, codon table replacement, engine architecture overhaul |
| **Minor** (`3.X.0`) | New rule, scan feature, optimization mode, new profile, new CLI flag |
| **Patch** (`3.1.X`) | Bug fix, metric correction, documentation update, dependency patch |

**Release checklist:**
1. Move `[Unreleased]` entries to `[X.Y.Z] — YYYY-MM-DD` in this file
2. Bump `version` in `pyproject.toml`
3. `git commit -m "chore: release vX.Y.Z"`
4. `git tag -a vX.Y.Z -m "Release vX.Y.Z"` → `git push && git push --tags`
5. Create GitHub Release from the tag

---

## [Unreleased]

---

## [3.1.0] — 2026-05-24

### Added

- **Custom restriction site removal** — synonymous substitution of user-specified restriction sites via API and web UI.
- **Rare codon run detection** (`scan_rare_codon_runs`) — detects consecutive rare codons (w < 0.3, default `min_run=3`) for ribosome stalling risk; included in full scan mode, excluded from fast mode.
- **Dinucleotide reduction pass** (`fix_dinucleotides`) — greedy CpG/TpA synonymous substitution integrated into the v2 optimization pipeline.
- **CAI-budgeted dinucleotide modes** — `aggressive` (no CAI check), `balanced` (CAI floor), `cai_preserving` (max CAI drop budget); `cai_before` / `cai_after` added to return dict.

### Fixed

- **Pipeline metric accuracy** — `candidate_metrics` CAI, GC%, and composite score are now re-measured from the final sequence after dinucleotide fix. Previously, pre-fix values were reported, understating CAI by ~0.04 on average (49/49 sequences now correctly report CAI ≥ 0.75).
- **CAI guard weight consistency** — `_calc_cai()` now uses golden-set reference weights (Sharp & Li 1987), consistent with `calculate_cai()`. Previously used working-table weights (~0.15 higher), which prevented the `balanced` guard from ever triggering.

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

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.1.0...HEAD
[3.1.0]: https://github.com/eijex/factorforge-cds/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/eijex/factorforge-cds/releases/tag/v3.0.0
