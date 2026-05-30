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

**Pre-release gate (before bumping):**
0. `git status --short` — working tree must be clean
0. `python -m ruff check .` — no lint errors
0. `python -m pytest tests/ -v --tb=short` — all tests pass
0. `python scripts/release.py X.Y.Z --dry-run` — verify 16 files

**Version bump & manual updates:**
1. Move `[Unreleased]` entries to `[X.Y.Z] — YYYY-MM-DD` in this file; update comparison links at bottom
2. Run `python scripts/release.py X.Y.Z` — updates all 16 version-bearing files + residual check
3. Add changelog entry to `web/index.html` (version panel HTML — manual; set new block to emerald/Current, demote previous to gray)
4. Add summary entry to `docs/changelog.md`

**Commit & CI gate (tag AFTER CI passes):**
5. `git commit -m "chore: release vX.Y.Z"` → `git push`
6. Wait for CI to pass: https://github.com/eijex/factorforge-cds/actions
7. `git tag -a vX.Y.Z -m "Release vX.Y.Z"` → `git push --tags`
8. GitHub Actions publishes to PyPI + Docker; creates GitHub Release → Zenodo DOI issued automatically

**Post-release verification:**
9. Verify PyPI: `pip install factorforge-cds==X.Y.Z && factorforge --help` (smoke test)
10. Verify Docker: `docker run ghcr.io/eijex/factorforge-cds:vX.Y.Z factorforge --help` (smoke test)
11. Confirm Zenodo DOI: https://zenodo.org/doi/10.5281/zenodo.20407331
12. **Bioconda** — update `recipes/meta.yaml` (version + SHA256 via `curl -s https://pypi.org/pypi/factorforge-cds/X.Y.Z/json | python -c "import sys,json; d=json.load(sys.stdin); [print(f['digests']['sha256']) for f in d['urls'] if f['packagetype']=='sdist']"`), push to fork branch `add-factorforge-cds`. Once PR is merged by Bioconda maintainers, autobump handles subsequent releases automatically.
13. **GitHub Issues** — close all issues completed in this release; if a full milestone is done, close the milestone too (`gh api repos/eijex/factorforge-cds/milestones/{N} --method PATCH --field state=closed`)

**Conditional checklists — apply only when relevant:**

<details>
<summary>New expression host added (e.g. BY-2, Arabidopsis)</summary>

- [ ] `src/factorforge/data/{host}_codons.json` — new codon table
- [ ] `src/factorforge/engines/profile/rules/` — host param plumbing
- [ ] `pyproject.toml` keywords — add new species names
- [ ] `CITATION.cff` title, abstract, keywords
- [ ] `web/index.html` `<title>` and `<meta description>`
- [ ] `README.md` tagline
- [ ] `docs/index.md` tagline + Supported Hosts table row
- [ ] `docs/how-it-works.md` — host-agnostic description
- [ ] `docs/cli.md` — `--host` option choices
- [ ] `docs/output.md` — CAI description
- [ ] `docs/validation.md` — GC% range note
- [ ] `docs/profiles.md` — Supported Hosts section
- [ ] `mkdocs.yml` `site_description`
- [ ] Documentation sync (internal) — update capability descriptions
- [ ] `_analysis/` — codon table validation analysis job (prerequisite)

</details>

<details>
<summary>New algorithm added (tAI, codon pair bias, 5' UTR MFE, etc.)</summary>

- [ ] `CHANGELOG.md` [Unreleased] entry
- [ ] Documentation sync (internal) — update differentiators
- [ ] `docs/how-it-works.md` — pipeline stage description
- [ ] `docs/profiles.md` — new scan/feature description

</details>

<details>
<summary>Bioconda PR merged (first time only)</summary>

- [ ] `docs/index.md` Access Options table — add conda install
- [ ] `docs/getting-started.md` — add `conda install -c bioconda factorforge-cds`
- [ ] `README.md` — add Bioconda badge + install instructions
- [ ] Documentation sync (internal) — add Bioconda distribution note

</details>

---

## Development History

FactorForge's public release history (v3.0+) builds on earlier internal implementation generations. Archived tracks are preserved under `archive/` for provenance and are not part of the installed package.

| Generation | Status | Description |
|-----------|--------|-------------|
| v1 — NBent_OptiCodon | Internal | Thesis-derived codon optimization baseline for *N. benthamiana* |
| v2 — Rule-Based Engine | Internal → Production | Deterministic constraint-aware engine; matured into `factorforge.engines.profile` |
| v3-alpha — ML Prototype | Archived | ML-based design attempt; performance insufficient; see `archive/v3-ml-prototype/` |
| v3.0+ — Current release | Public | Open-source release of the v2 engine; development continues here |
| v3.7+ — ML Engine | Planned | ML-based design as `--engine ml`; added once sufficient wet-lab data is available |

---

## [Unreleased]

---

## [3.1.6] — 2026-05-30

### Added

- **SynCodonLM scoring dimension** — optional 5th composite score component (`w_syncodonlm`, default `0.0`). Integrates Boehringer-Ingelheim's BERT-based codon language model ([SynCodonLM, NAR 2025](https://github.com/Boehringer-Ingelheim/SynCodonLM); HuggingFace: `jheuschkel/SynCodonLM-V2`). Graceful fallback (score 0.5, WARNING) when `transformers` is not installed. No change to existing scoring behavior.
- **`ml_enhanced` scoring profile** — `w_cai=0.35, w_gc=0.25, w_mfe=0.15, w_syncodonlm=0.25`. Opt-in; existing four profiles unchanged.
- **`[ml]` optional dependency group** — `pip install factorforge-cds[ml]` installs `transformers>=4.40` and `torch>=2.0` for SynCodonLM inference.
- **`scoring_ml.py`** — `SynCodonLMScorer` class with lazy model loading; `calculate_syncodonlm_score(sequence, organism)`.
- **Profile comparison mode** — `factorforge optimize input.fasta --engine profile --compare-profiles balanced,high_cai,gc_target` outputs a side-by-side CAI / GC% / score table. First profile result saved to `--output` when specified. `POST /api/optimize/compare` endpoint added with same functionality via JSON API.
- **Tutorial: GFP N. benthamiana** — end-to-end worked example at `docs/tutorials/gfp-nbenthamiana.md`. Covers CLI, Python API, profile comparison, and MoClo assembly preparation.
- **Batch optimization API** — `POST /api/optimize/batch` accepts up to 20 sequences in a single request. Returns per-sequence CAI, GC%, score, and optimized CDS. Auto-generates IDs (`seq_1`, `seq_2`, ...) when omitted. CLI multi-FASTA was already supported.
- **Tobacco BY-2 host support (experimental)** — `--host by2` CLI flag and `"host": "by2"` API field optimize for *N. tabacum* BY-2 suspension culture cells using a Kazusa-derived codon table (1,534 CDS, species 4097). Default host remains `nbenthamiana`. CAI difference between hosts is < 0.05. **Experimental:** uses *N. tabacum* codon usage as proxy; not wet-lab validated for BY-2 expression performance.
- **Structure prediction links** — AlphaFold DB and ESM Atlas fold links appear in the result panel after optimization. No API calls — links open external services with the input sequence.

---

## [3.1.5] — 2026-05-28

### Fixed
- **Data file packaging** — `pyproject.toml`에 `package-data` 누락으로 `pip install` 후 `nbenthamiana_codons.json` 경로 오류 발생하던 버그 수정. PyPI wheel에 JSON 데이터 파일이 포함되지 않았던 문제.
- **Path resolution** — `domesticator.py`, `pipeline.py`, `analysis/metrics.py`의 `parents[N]` 하드코딩을 `get_data_path()` 사용으로 통일. 설치 환경에서 경로가 `Lib/data/`로 잘못 해석되던 문제 수정.
- **CI matrix** — Python 3.13 + Windows (`windows-latest`) 테스트 추가 (3.10/3.11/3.12/3.13 × ubuntu/windows 8개 조합).

---

## [3.1.4] — 2026-05-27

### Added

- **CITATION.cff** — GitHub "Cite this repository" button; updated at every version bump.
- **SECURITY.md** — vulnerability reporting policy (GitHub Private Vulnerability Reporting + email).
- **ROADMAP.md** — public development direction: v3.2 planned features, ML Research Track (v4), wet-lab scope.
- **bump_version.py** — automates version string updates across 14 files (`python scripts/release.py X.Y.Z`).
- **Development history narrative** — README, CHANGELOG, archive READMEs, and docs/changelog.md now document the v1→v2→v3-alpha→v3.x→v4 version lineage.

### Changed

- **Public engine naming** — promoted the stable rule/profile engine to `factorforge.engines.profile` and CLI `--engine profile`.
- **Package metadata** — switched to SPDX license metadata and removed the stale archived v1 CLI runtime path from the installed package surface.
- **Benchmark numbers corrected** — docs/index.md updated with v3.1.3 measurements (N=3,876 SGN CDS, balanced profile): CAI 0.76, GC% 59.77%, target range 55–65%. Previous values were measured with pre-v3.1.0 parameters.

### Fixed

- **Stale test version** — `tests/test_schemas/test_design_package.py` `product_version` bumped from 3.1.1 → 3.1.3.

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

- **Profile Rule-based Optimizer** — CAI-maximizing baseline with profiles: `balanced`, `high_cai`, `gc_target`, `viral_delivery` (`src/factorforge/engines/profile/`)
- **DP Feasibility Engine** — constraint-based dynamic programming; computes achievable CAI/GC ranges under synonymous codon constraints
- **Metrics Engine** — CAI, GC (global/local/first-region), amino acid identity, internal stop count, homopolymer/repeat/forbidden motif/invalid codon detection
- **Sequence Validator** — structured hard-fail contract (AA identity, stops, invalid codons, GC, length) returning machine-readable dict (`src/factorforge/utils/validation.py`)
- **SGN data pipeline** — fetch, CAI-filter, and split *N. benthamiana* CDS from Sol Genomics Network v2.6.1
- **Benchmark Panel** — 5 synthetic proteins in `tests/fixtures/benchmark_proteins.py`
- **Diagnostic & benchmark scripts** — reproducible sequence metric checks and feasibility benchmarks
- **CLI** — `factorforge optimize` / `factorforge list-engines` with FASTA and GenBank output
- **Docs** — public installation, CLI, profile, output, validation, and changelog pages

### Tests

- 333 passing at release time

---

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.1.6...HEAD
[3.1.6]: https://github.com/eijex/factorforge-cds/compare/v3.1.5...v3.1.6
[3.1.5]: https://github.com/eijex/factorforge-cds/compare/v3.1.4...v3.1.5
[3.1.4]: https://github.com/eijex/factorforge-cds/compare/v3.1.3...v3.1.4
[3.1.3]: https://github.com/eijex/factorforge-cds/compare/v3.1.2...v3.1.3
[3.1.2]: https://github.com/eijex/factorforge-cds/compare/v3.1.1...v3.1.2
[3.1.1]: https://github.com/eijex/factorforge-cds/compare/v3.1.0...v3.1.1
[3.1.0]: https://github.com/eijex/factorforge-cds/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/eijex/factorforge-cds/releases/tag/v3.0.0
