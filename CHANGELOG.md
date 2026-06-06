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
0. Verify `Dockerfile` references match current file locations (especially after any script/file moves)
0. `python scripts/release.py X.Y.Z --dry-run` — verify 16 files
0. `python scripts/release.py X.Y.Z --dry-run` — verify cross-repo docs + MCP tool strings

**Version bump & manual updates:**
1. Move `[Unreleased]` entries to `[X.Y.Z] — YYYY-MM-DD` in this file; update comparison links at bottom
2. Run `python scripts/release.py X.Y.Z` — updates all version-bearing files + cross-repo docs + MCP tool strings + residual check
3. Add changelog entry to `web/index.html` (version panel HTML — manual; set new block to emerald/Current, demote previous to gray)
4. Add summary entry to `docs/changelog.md`
5. Review internal workspace docs (roadmap, job catalog, memory) for any updates needed

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
13. **GitHub house maintenance:**
    - [ ] Close issues completed in this release; close milestone if all items are done (`gh api repos/eijex/factorforge-cds/milestones/{N} --method PATCH --field state=closed`)
    - [ ] Close stale auto-generated CI failure issues (`gh issue list --label ci-failure` → close any that are no longer active)
    - [ ] Scan open issues for policy violations (paper/JOSS references, unpublished roadmap items) → close any found
    - [ ] ROADMAP.md version numbers up to date (check if any milestone version labels changed)
    - [ ] Issue template (`wet_lab_result.yml`) fields match Google Form? → sync if diverged
    - [ ] VALIDATION.md current? → update GC% ranges, host list, "How to Contribute" fields if needed
14. **eijex-mcp sync** — check the eijex-mcp repo for required updates:
    - [ ] Version string in `factorforge_cds_optimize` description → bump to new version
    - [ ] New profiles added? → add to `profile` enum in route.ts and mcp-tools.ts
    - [ ] New API endpoints added? → expose as new MCP tools
    - [ ] New `--host` options? → update tool description
    - [ ] `CHANGELOG.md` — add entry for any tool additions or interface changes
    - [ ] `README.md` tool table still accurate?
    - [ ] `skills/bio-research.md` workflow steps still accurate?
    - [ ] `src/app/page.tsx` section labels/descriptions + safety disclaimer still accurate?
    - [ ] `src/app/_lib/mcp-tools.ts` tool descriptions, tags, keyFeatures still accurate?
    - [ ] `.claude-plugin/plugin.json` description still accurate?
    - [ ] `package.json` version — bump eijex-mcp if interface changed
    - [ ] `CHANGELOG.md` — add entry for FactorForge version bump and any breaking changes affecting MCP tools
15. **eijex-web** — update version display (when repo is active):
    - [ ] `src/app/components/StatsBar.tsx` — `latest release` value
    - [ ] `src/app/components/Footer.tsx` — `factorforge-cds vX.Y.Z` string
16. **Bioconda PR** — update PR title to new version: `gh pr edit 65834 --repo bioconda/bioconda-recipes --title "Add factorforge-cds X.Y.Z"` + add comment with new SHA256
16. **Wet-lab survey** — check if Google Form host/profile fields need updating when a new host or profile is added
17. **www.eijex.com sync (eijex-web)** — **conditional**, only when the release changes something the company site shows. The site links to FactorForge (factorforge.eijex.com, `pip install`, Docker, Bioconda) without pinning a version, so routine version bumps need **no** update. Update `C:\Work\eijex\eijex-web` only when:
    - [ ] New capability/claim worth showcasing on the product/landing copy
    - [ ] New install channel (e.g. Bioconda first becomes available) → add to install links
    - [ ] A hard-coded FactorForge version string was introduced anywhere on the site
    - [ ] Vercel auto-deploys eijex-web on push to its repo — no tag needed

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
- [ ] Internal codon table validation analysis (prerequisite)
- [ ] **Google Form** — add new host to "Host organism" field options

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
| v4.0 — ML Engine | Planned | ML-based design as `--engine ml`; added once sufficient wet-lab data is available |

---

## [Unreleased]

### Added
- **MFE metadata fields** — Design Package and API response now include `mfe_used` (bool), `mfe_status` (`computed` / `not_computed`), and `mfe_warning` (string when ViennaRNA unavailable). `score_components` added to expose per-term weights used in composite score calculation.

### Fixed
- **Domestication Silence Fail** — `pipeline.py` now raises `ValueError` when restriction-site domestication fails (previously returned the undomesticated sequence silently as success).
- **Pipeline Output Validator** — `validate_cds_output()` is now called in `pipeline.py` before final sequence return, catching AA identity violations and internal stops at the pipeline level (previously only enforced at the API layer since job 067).
- **MFE not-computed value** — `mfe_kcal_mol` is now `null` (not `0.0`) when ViennaRNA is unavailable (e.g. production Vercel). Composite score is unchanged; this corrects misleading metadata only.

### Documentation
- **Stale constant corrections** — 5 doc/comment locations corrected to match live code: `target_cai` 0.92 → 0.82; GC feasibility range 41–44% → 55–65%; 5′ ramp docstring bottom 50% → 25%; `scan_gc_extremes` local-guard behaviour documented; feasibility fallback path description corrected.

---

## [3.1.9] — 2026-06-04

### Documentation
- **Internal housekeeping** — project tracking references updated. No engine changes.

---

## [3.1.8] — 2026-06-01

### Breaking Changes
- **`gc_target` profile default changed** — calling `gc_target` without an explicit `target_gc` now produces sequences targeting ~60% GC (host midpoint) instead of the previous 42.5%. If you relied on the 42.5% default, pass `target_gc=42.5` explicitly to preserve the old behavior.

### Changed
- **`gc_target` profile default** — now targets the host-profile GC midpoint (60% for *N. benthamiana*) when `target_gc` is not supplied, instead of the legacy hardcoded 42.5%. To target a lower GC, pass `target_gc` explicitly. Output sequences from `gc_target` without an explicit target will differ.
- **GC scoring** — `calculate_composite_score` now scores GC using a band function (`gc_band_score`): full score inside `[gc_min, gc_max]`, linear decay outside over `gc_decay_width` (default 20 pp). Replaces the previous `1 - |GC - GC_opt|/50` proximity formula, which under-discriminated GC quality.
- **`assembly_friendly` scoring weights** — changed from balanced-identical `(0.5, 0.3, 0.2)` to `(0.3, 0.4, 0.3)` (lower CAI pressure, higher GC/MFE weight) to align scoring with its Type IIS site-avoidance translation strategy.
- **`feasibility.py` defaults** — `target_cai` 0.92 → 0.82 (achievable; aligns with industry >0.8 practice); `target_gc` 41–44% → 55–65%; fallback GC ranges realigned to the 55–65% output distribution.

### Fixed
- **Homopolymer thresholds documented** — expression-stability (≥6 nt) and synthesis/manufacturing (≥8 nt) scans now use named constants (`HOMOPOLYMER_EXPRESSION_WARN_NT`, `HOMOPOLYMER_SYNTHESIS_WARN_NT`) and emit `context`/`threshold_nt` metadata so the two intentionally different thresholds are no longer mistaken for a bug.
- **Misleading docs removed** — `gc_target` no longer described as "42.5% (N. benthamiana optimal)"; 42.5% was a legacy assumption inconsistent with the 55–65% codon-table output.
- **CLI docs corrected** — `docs/cli.md` `--gc-min`/`--gc-max` defaults fixed from 40/55 to the actual 55/65.

### Documentation
- **`docs/profiles.md`** — added missing `assembly_friendly` profile; corrected `gc_target` description.
- **`docs/tutorials/gfp-nbenthamiana.md`** — regenerated profile-comparison metrics under the new GC scoring and `gc_target` default.

---

## [3.1.7] — 2026-05-31

### Added
- **Web UI host selector** — expression host toggle (N. benthamiana / BY-2 Experimental) in the input panel. BY-2 selection disables Feasibility Best objective and shows experimental warning. Result panel displays active host profile.
- **E2E smoke tests** — 5 Playwright smoke tests covering UI load, protein input, BY-2 host toggle, Feasibility Best guard, and result rendering. Runs automatically after each deployment via `e2e.yml`.

### Documentation
- **Eijex MCP access** — added Eijex MCP as access option in `README.md` and `docs/index.md`
- **API endpoints** — added `POST /api/optimize`, `/compare`, `/batch` endpoints section to `docs/cli.md`
- **MCP getting started** — added Eijex MCP connection guide to `docs/getting-started.md`
- **ml_enhanced profile** — `docs/profiles.md`에 ml_enhanced 프로파일 문서화
- **AGENTS.md** — 새 API 엔드포인트 추가 시 eijex-mcp 동기화 항목 명시

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
- **Tobacco BY-2 host support (experimental)** — `--host by2` CLI flag and `"host": "by2"` API field optimize for *N. tabacum* BY-2 suspension culture cells using a Kazusa-derived codon table (1,534 CDS, species 4097). Default host remains `nbenthamiana`. **Experimental:** uses *N. tabacum* codon usage as proxy; not wet-lab validated for BY-2 expression performance.
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

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.1.9...HEAD
[3.1.9]: https://github.com/eijex/factorforge-cds/compare/v3.1.8...v3.1.9
[3.1.8]: https://github.com/eijex/factorforge-cds/compare/v3.1.7...v3.1.8
[3.1.7]: https://github.com/eijex/factorforge-cds/compare/v3.1.6...v3.1.7
[3.1.6]: https://github.com/eijex/factorforge-cds/compare/v3.1.5...v3.1.6
[3.1.5]: https://github.com/eijex/factorforge-cds/compare/v3.1.4...v3.1.5
[3.1.4]: https://github.com/eijex/factorforge-cds/compare/v3.1.3...v3.1.4
[3.1.3]: https://github.com/eijex/factorforge-cds/compare/v3.1.2...v3.1.3
[3.1.2]: https://github.com/eijex/factorforge-cds/compare/v3.1.1...v3.1.2
[3.1.1]: https://github.com/eijex/factorforge-cds/compare/v3.1.0...v3.1.1
[3.1.0]: https://github.com/eijex/factorforge-cds/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/eijex/factorforge-cds/releases/tag/v3.0.0
