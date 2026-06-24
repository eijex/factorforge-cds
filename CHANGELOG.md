# Changelog

All notable changes to FactorForge are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Release Policy

| Bump | When to use |
|------|-------------|
| **Major** (`X.0.0`) | Breaking API change, codon table replacement, engine architecture overhaul |
| **Minor** (`3.X.0`) | New rule, scan feature, optimization mode, new profile, new CLI flag |
| **Patch** (`3.1.X`) | Bug fix, metric correction, documentation update, dependency patch |

Before each release, maintainers run lint/tests, update version-bearing files,
verify package and Docker smoke tests, and check public documentation for
version drift, unsupported claims, sensitive-data guidance, and stale examples.

---

## [Unreleased]

### Added

- Runtime validation registry (`factorforge.validation_registry`) listing all 17 validation checks with per-execution-path enforcement metadata, and a canonical validation report builder (`factorforge.validation_report`) that scans the final returned CDS and reports results keyed by `check_id`. The web "Sequence Checks" panel and the `/api/optimize` response now expose all 9 advisory scanners plus restriction-site and MoClo-overhang results (previously only 3 of 9+ checks were visible). `validation.moclo`/`polya`/`gc` legacy fields are unchanged.
- `reproducibility/benchmark_v0.5.1/scripts/figures/make_cai_vs_multiconstraint_figure.py` regenerates the "CAI-only baseline versus constraint-aware design" scatter figure (mean CAI vs. corrected multi-constraint pass rate per method) from `benchmark_summary.frozen.json`. No generation script for this figure existed anywhere in the repository or its Git history; this is a reconstruction verified against the figure's existing per-method values, not a recovery of the original code.

### Docs

- Clarified that BY-2 (N. tabacum) experimental host support applies to
  the `balanced`/`gc_target`/`assembly_friendly` rule-based profiles only;
  `high_cai` and `feasibility_best` remain N. benthamiana-only by design.
- `docs/how-it-works.md`'s "Design Objectives" table listed `gc_target`/
  `high_cai` as if they were working DP `--objective` values; passing either
  to `--objective` always raised (the DP engine only implements
  `feasibility_best`). Split the table into the DP engine's one real
  objective and the profile engine's `--profile` values, with an explicit
  note that the latter are not DP objectives.
- `docs/profiles.md`'s "Stable Profiles" section claimed all four profiles
  are "fully supported ... via CLI, Python API, web app, and MCP" without
  disclosing that `high_cai` is N. benthamiana-only and rejected for other
  hosts at the CLI/REST/web layers. Added that qualifier.

### Fixed

- `/api/optimize/compare` and `/api/optimize/batch` accepted a `host` or
  `host_profile` field in the request body but never read it, silently
  ignoring the caller's host intent and returning HTTP 200 with
  default-host output. Both endpoints now reject any request containing
  either field with HTTP 400 (`HOST_NOT_SUPPORTED_ON_ENDPOINT`) instead of
  silently dropping it.
- CLI `--objective` declared `gc_target` and `high_cai` as valid DP
  objective values via `click.Choice`, but the DP engine only ever
  implements `feasibility_best` — passing either always raised
  `DP engine currently supports --objective feasibility_best.`. Removed
  both from the `click.Choice` so the CLI rejects them at the argument-parsing
  stage with a clear "invalid choice" error instead of failing inside the
  engine.
- A direct library call (`RuleBasedOptimizer().optimize(profile="high_cai",
  host=<non-default>)`) silently substituted N. benthamiana golden-set
  output with no indication that the requested host was ignored (the
  existing host/strategy compatibility guard only covers the REST/CLI/web
  surfaces). Added a `logger.warning` at this call site; the host-invariant
  output itself is unchanged (this is the documented design boundary, not a
  bug) and the existing benchmark
  `codon_table_path` injection path is unaffected.
- Web UI "View Predicted Structure → ESM Atlas" link used a URL scheme
  (`esmatlas.com/explore?tab=fold&sequence=`) that ESM Atlas's router never
  branches on, so it always landed on the generic explore screen instead of
  a fold view regardless of sequence. Switched to the route ESM Atlas's own
  navigation actually uses (`esmatlas.com/resources?action=fold&sequence=`),
  confirmed against ESM Atlas's production JS bundle and verified live —
  it now lands on the correct "Fold Sequence" tool (ESM Atlas's UI does not
  support pre-filling the sequence box from a URL parameter, so users still
  paste it manually after arriving).
- Removed internal task-tracking references ("Job NNN") from the public changelog (`CHANGELOG.md`, `docs/changelog.md`, `web/index.html`) and from public docs (`docs/strategy/eijex-tool-layer-classification.md`, `docs/validation/RELEASE_GATE.md`, `reproducibility/benchmark_v0.5.1/README.md`). `scripts/audit_public_surface.py`'s `internal_reference` pattern now also matches `Job \d+` so future releases catch this automatically.
- Removed the same internal task-tracking references ("Job NNN", "analysis NNN") from source-level comments, docstrings, and registry provenance fields (`src/`, `tests/`), and deleted `benchmarks/scripts/resume_job130_rerun.sh`, a one-off maintainer rerun script that hardcoded a local developer machine path. Empirical citations (CAI/GC benchmark numbers) were kept; only the internal ID labels were removed.
- `/api/optimize` now rejects explicit `objective=feasibility_best` or
  `profile=high_cai` requests combined with a non-default `host` (HTTP
  400) instead of silently returning N. benthamiana-table output. The
  implicit case (host-only request, no explicit strategy) now discloses
  `requested_strategy`/`resolved_strategy`/`resolution_reason` and
  resolves to `balanced` (previously silently resolved toward `high_cai`,
  which is also N. benthamiana-only). CLI and web UI now block `high_cai`
  for non-default hosts (`--compare-profiles` included), matching the
  existing `feasibility_best` guard; the web UI's auto-selected fallback
  for BY-2 changes from `high_cai` to `gc_target`.
- Release provenance hashing now computed from the committed git blob (`git show HEAD:<path>`) instead of local working-tree bytes, fixing CRLF/LF drift on Windows that could silently produce incorrect SHA-256 values in `reproducibility/benchmark_v0.5.1/MANIFEST.json` and `tests/test_docs_consistency.py`.
- Public-surface DOI references (`README.md`, `docs/index.md`, `AGENTS.md`) switched from version-pinned Zenodo DOIs to the concept DOI, which always resolves to the latest release, so future releases no longer require manually updating these files.
- `AGENTS.md` no longer states a hardcoded "16 version-bearing files" count, which had gone stale; points at `scripts/release.py`'s `build_targets()` as the source of truth instead.
- `docs/rule-engine-roadmap.md`, `docs/validation.md`, `docs/how-it-works.md`, and `docs/factorforge-architecture.md` enumerated only 5-8 of the 9 default advisory RuleEngine scanners and omitted the MoClo overhang assembly-review check; all four now list the complete set. `rule-engine-roadmap.md` additionally mis-stated "Repeat patterns" as "Planned / Not yet implemented" (it is implemented and runs by default) and described an unused legacy GC-window calculation instead of the active `scan_gc_extremes` thresholds; both corrected. No runtime code changed.
- Web app "Sequence Checks" badge labeled "MoClo Overhang Check" actually reported a Type IIS restriction-site scan result, not MoClo overhang validity (the real MoClo overhang check lives in the opt-in construct-builder path and was never called here); label and changelog text corrected to "Restriction Site Check (Type IIS)". `validation.moclo` JSON field name kept unchanged for frontend compatibility.

### Added

- `scripts/audit_public_surface.py` and a CHANGELOG duplicate-`[Unreleased]`-header check now run on every push/PR in CI, instead of only when `release.py --auto --audit-script` is remembered.
- `scripts/regen_manifest.py` regenerates `reproducibility/benchmark_v0.5.1/MANIFEST.json` input hashes from committed git-blob content on demand.

### Changed

- Repositioned FactorForge as a claim-bounded pre-synthesis review harness across `README.md`, `ROADMAP.md`, `docs/index.md`, the new `docs/roadmap.md`, and `web/index.html`, separating the research-software journey from product roadmap themes without adding any new guarantees (expression, glycosylation, folding, yield, synthesis acceptance, regulatory approval).

---

## [3.2.2] — 2026-06-18

### Fixed

- **`multi_constraint_pass` definition corrected (scoring_contract v1.1)**: `benchmarks/scoring.py`
  `score_cds()` now defines `multi_constraint_pass = biological_pass AND assembly_pass AND gc_in_target_range`.
  The previous definition (`biological_pass AND assembly_pass`) omitted GC target compliance, producing
  inflated L3/L4 ablation values (89.0%/88.6%) that were mathematically inconsistent with their GC in range
  rates (3.7%/5.8%). Corrected values: L3=3.5%, L4=5.6%.
  All benchmark artifacts (`benchmark_summary.json`, `ablation_summary.json`, `benchmark_v0.5.1` data/figures)
  regenerated from full rerun (seed=320, N=49,257). Zenodo `benchmark_results.csv` v2 (DOI: 10.5281/zenodo.20676276) supersedes v1.
  A `canonical_multi_constraint_pass()` helper added for recomputing from primitive columns in historical CSVs.
- Benchmark: source-profile codon-table injection now flows into both the design and scoring paths (`design_table_sha256 == score_table_sha256` verified per run), fixing a prior gap where design always used the default table regardless of an injected profile.

### Added

- Data: added three genome-annotated *N. benthamiana* codon-usage profiles (SGN QLD183 v103 CDS-derived; SGN NbeV1.1 all-CDS-derived; SGN NbeV1.1 high-confidence-CDS-derived) built via `scripts/build_codon_profile.py` under `strict_nuclear_cds_v1` filtering, alongside the existing packaged reference profile.

### Changed

- Web: Host System cards in `web/index.html` are now rendered dynamically from `GET /api/optimize` (`supported_hosts` + new `host_metadata` field) instead of being hardcoded, removing a 3-way duplication between the HTML, `web/js/app.js`, and the API.

### Docs

- Aligned public wet-lab validation contribution language with manual-review, public-safe submission rules.
- Clarified that public GitHub Issues must not contain raw sequences, confidential construct details, internal batch IDs, patient data, private contact information, exact process parameters, or confidential partner/customer data.
- Aligned public README, docs, web, citation, packaging, roadmap, and benchmark wording with the in-silico CDS design claim boundary.
- Removed maintainer-local file paths and internal repo references from `docs/release-checklist-template.md`, replacing them with generic placeholders.

---

## [3.2.1] — 2026-06-16

### Added
- Protein risk annotation layer for CDS sequences
  - Transmembrane helix prediction (Kyte-Doolittle, window=19, threshold=1.6)
  - Signal peptide heuristic (N-terminal 30 aa scan)
  - Risk classification: HIGH / MEDIUM / LOW / UNKNOWN

### Fixed
- Correct CAI provenance annotation in benchmark output
- Correct Type IIS restriction site warning status
- Pin JSON files to LF line endings for Windows reproducibility
- Add pandas to dev test dependencies
- Manifest SHA-256 reproducibility drift on Windows (JSON/EOL normalization)
- wet-lab result GitHub template: add protein_class options (Reporter / Antigen / Cytokine) and validation consent checkbox

### Changed
- Standardize wet-lab submission link labels across README, docs, and web app

### Docs
- Add public-safe FactorForge agent guidance
- Expand release checklist with public surface audit steps

### Chore
- Bump actions/configure-pages, codecov/codecov-action, actions/deploy-pages, actions/setup-python, softprops/action-gh-release (Dependabot)

---

## [3.2.0] — 2026-06-11

### Added
- **MFE metadata fields** — Design Package and API response now include `mfe_used` (bool), `mfe_status` (`computed` / `not_computed`), and `mfe_warning` (string when ViennaRNA unavailable). `score_components` added to expose per-term weights used in composite score calculation.
- **Design Package schema v1.0.0** — Formal IUPAC/FASTA I/O contracts and MFE null invariant established (090).
- **Registry production constants export** — `DEFAULT_CAI_TARGET`, `DEFAULT_GC_LOW`, `DEFAULT_GC_HIGH` importable as public production constants (091).
- **Benchmark seed injection** — `--seed` flag for deterministic reruns; `most_frequent_codon` tie-breaking deduplication (099).
- **Codon table provenance disclosure** — `codon_table_manifest.json` with sha256 pin, `build_path_status: incomplete`, and known limitations for `nbenthamiana_codons.json` (097).

### Fixed
- **Domestication Silence Fail** — `pipeline.py` now raises `ValueError` when restriction-site domestication fails (previously returned the undomesticated sequence silently as success).
- **Pipeline Output Validator** — `validate_cds_output()` is now called in `pipeline.py` before final sequence return, catching AA identity violations and internal stops at the pipeline level.
- **MFE not-computed value** — `mfe_kcal_mol` is now `null` (not `0.0`) when ViennaRNA is unavailable. Composite score is unchanged; this corrects misleading metadata only.
- **Input validator** — IUPAC ambiguous DNA/AA sequence misclassification corrected (098).

### Documentation
- **Stale constant corrections** — 5 doc/comment locations corrected to match live code.
- **Claim wording alignment** — Public-facing API and CLI output wording unified; no expression-level or yield improvement claims (092).
- **Formal benchmark** — N. benthamiana SGN CDS (N=49,257, seed=320). All metrics are in-silico; no wet-lab validation claimed.

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
- **Profile documentation** — public profile table clarified supported and internal profile boundaries

---

## [3.1.6] — 2026-05-30

### Added

- **Experimental scoring hooks** — optional internal scoring hooks were added without changing existing public profile behavior.
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
- **ROADMAP.md** — public development direction, validation scope, and planned host/profile work.
- **bump_version.py** — automates version string updates across 14 files (`python scripts/release.py X.Y.Z`).
- **Public history cleanup** — archive references kept outside the primary user docs.

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

- **Wet-lab feedback modal** — "Submit Wet-lab Result" button opens the public-safe GitHub feedback template.
- **JSON Copy button** — one-click copy of the full optimization JSON output in the terminal panel.
- **Submit button tooltip** — info icon explains what wet-lab feedback is used for.

### Changed

- **Design Objective order** — reordered to match recommended wet-lab testing sequence: Feasibility Best → 5' Ramp → High CAI → GC Target → Assembly Friendly → Viral Delivery.
- **Wet-lab feedback fields** — GitHub Issue Template updated for structured public-safe feedback.
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

[Unreleased]: https://github.com/eijex/factorforge-cds/compare/v3.2.2...HEAD
[3.2.2]: https://github.com/eijex/factorforge-cds/compare/v3.2.1...v3.2.2
[3.2.1]: https://github.com/eijex/factorforge-cds/compare/v3.2.0...v3.2.1
[3.2.0]: https://github.com/eijex/factorforge-cds/compare/v3.1.9...v3.2.0
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
