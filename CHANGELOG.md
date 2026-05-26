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
1. Move `[Unreleased]` entries to `[X.Y.Z] — YYYY-MM-DD` in this file; update comparison links at bottom
2. Bump `version` in `pyproject.toml`
3. Update version strings: `src/factorforge/__init__.py`, `src/factorforge/engines/profile/__init__.py`, `src/factorforge/engines/v2/__init__.py`, `src/factorforge/engines/v3/__init__.py`, `src/factorforge/engines/__init__.py`, `src/factorforge/engines/profile/optimizer.py`, `src/factorforge/engines/v3/pipeline.py`, `api/optimize.py` (comment + ENGINE_VERSIONS), `web/index.html` (button + changelog entry), `web/js/app.js`
4. Update tests: `tests/api/test_optimize_contract.py`, `tests/engines/profile/test_profile_engine_alias.py`, `tests/engines/v2/test_cli_optimize.py`
5. Update `README.md` version badge and citation
6. Update `docs/index.md` version badge
7. Update `docs/changelog.md` — add new version entry (mirrors CHANGELOG.md, summarized)
6. `git commit -m "chore: release vX.Y.Z"`
7. `git tag -a vX.Y.Z -m "Release vX.Y.Z"` → `git push && git push --tags`
8. Create GitHub Release from the tag

---

## [Unreleased]

### Changed

- **Profile engine namespace** — promoted the stable rule/profile engine to `factorforge.engines.profile` and CLI `--engine profile`; retained `factorforge.engines.v2` and `--engine v2` as compatibility aliases.

---

## [3.1.3] — 2026-05-26

### Fixed

- **Disabled profile cards** — removed `pointer-events-none` so tooltips are accessible on hover; added "Pending wet-lab validation before enabling" notice to 5' Ramp and Viral Delivery cards.
- **Viral Delivery tooltip** — corrected citation reference; updated to Peccoud et al. 2024 (PMC11718241).
- **Analytics notice** — "Sequence content is not intentionally stored" → "Submitted sequences are not logged or stored" for clarity.

---

## [3.1.2] — 2026-05-26

### Fixed

- **viral_delivery scoring** — corrected citation reference; adjusted weights: `w_gc` 0.25→0.35, `w_mfe` 0.40→0.30 per PMC11718241 (Peccoud et al. 2024).
- **5' Ramp deoptimization** — N-terminal codon deoptimization strength reduced from bottom 50% to bottom 25% of frequency-sorted codons; aligns with PMC11718241 tAI_ramp 0.8–1.2 optimal range.

### Changed

- **Changelog label** — "ML Research Track" renamed to "Research Track".

---

## [3.1.1] — 2026-05-24

### Added

- **Wet-lab feedback modal** — "Submit Wet-lab Result" button opens an embedded Google Form iframe with version and optimization profile pre-filled from the current result.
- **JSON Copy button** — one-click copy of the full optimization JSON output in the terminal panel.
- **Submit button tooltip** — info icon explains what wet-lab feedback is used for.

### Changed

- **Design Objective order** — reordered to match recommended wet-lab testing sequence: Feasibility Best → 5' Ramp → High CAI → GC Target → Assembly Friendly → Viral Delivery.
- **Wet-lab feedback fields** — GitHub Issue Template updated with promoter, subcellular targeting, harvest timepoint, and native control fields.
- **VALIDATION.md / docs/validation.md** — "Include" section updated with new experimental metadata fields.

### Fixed

- **Vercel deployment** — `/api/optimize` was returning 404 due to Root Directory being set to `web/` in Vercel project settings; resolved by clearing Root Directory to repo root.

---

## [3.1.0] — 2026-05-24

### Added

- **Custom restriction site removal** — synonymous substitution of user-specified restriction sites via API and web UI.
- **Rare codon run detection** (`scan_rare_codon_runs`) — detects consecutive rare codons (w < 0.3, default `min_run=3`) for ribosome stalling risk; included in full scan mode, excluded from fast mode.
- **Dinucleotide reduction pass** (`fix_dinucleotides`) — greedy CpG/TpA synonymous substitution integrated into the profile optimization pipeline.
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

- **Profile Rule-based Optimizer** — CAI-maximizing baseline with profiles: `balanced`, `high_cai`, `gc_target`, `viral_delivery` (`src/factorforge/engines/profile/`; `src/factorforge/engines/v2/` remains a compatibility alias)
- **DP Feasibility Engine** — constraint-based dynamic programming; computes achievable CAI/GC ranges under synonymous codon constraints (`src/factorforge/ml/feasibility.py`)
- **Metrics Engine** — CAI, GC (global/local/first-region), amino acid identity, internal stop count, homopolymer/repeat/forbidden motif/invalid codon detection (`src/factorforge/ml/metrics.py`)
- **Sequence Validator** — structured hard-fail contract (AA identity, stops, invalid codons, GC, length) returning machine-readable dict (`src/factorforge/utils/validation.py`)
- **Synonym Mask** — per-position boolean mask restricting decoding to valid synonymous codons (`src/factorforge/engines/v3/synonym_mask.py`)
- **Constrained Decoder** — synonym-constrained autoregressive decoding + post-validation + profile engine fallback on failure (`src/factorforge/engines/v3/inference/constrained_decoder.py`)
- **Profile Engine Adapter** — formal baseline wrapper, returns structured dict with metrics and validator result (`src/factorforge/engines/v3/inference/v2_adapter.py`)
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

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.1.3...HEAD
[3.1.3]: https://github.com/eijex/factorforge-cds/compare/v3.1.2...v3.1.3
[3.1.2]: https://github.com/eijex/factorforge-cds/compare/v3.1.1...v3.1.2
[3.1.1]: https://github.com/eijex/factorforge-cds/compare/v3.1.0...v3.1.1
[3.1.0]: https://github.com/eijex/factorforge-cds/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/eijex/factorforge-cds/releases/tag/v3.0.0
