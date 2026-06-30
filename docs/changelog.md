# Changelog

Full changelog: [CHANGELOG.md on GitHub](https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md)

FactorForge v3.0+ is the current public release line.

## Unreleased

### Changed
- Activated NbeV1.1 high-confidence CDS-derived codon usage as the current *N. benthamiana* software default with a 40-47% GC policy band. Public reference overrides remain disabled; the legacy Kazusa/SGN composite is retained as a historical comparator. This is an in-silico software-default change only, not wet-lab validation, expression/yield prediction, or a biological-superiority claim.

## v3.2.8 — 2026-06-29

### Reverted
- Default *N. benthamiana* GC reference band / codon-reference table reverted to legacy (55-65% GC), undoing a default-codon-reference promotion that shipped in v3.2.7. Newer candidate/reference-comparator assets remain packaged for provenance and research comparison, but are not public production defaults or selectable UI options pending an MFE re-sensitivity + 2x2 factorial recheck

## v3.2.7 — 2026-06-29

### Changed
- Recalibrated the MFE-score normalization clamp from `[-0.5, 0.0]` to `[-0.40, -0.15]` based on a measured MFE distribution (n=135) showing the wider clamp was range-compressing scoring discriminability

### Fixed
- `calculate_composite_score()` now raises on an unrecognized profile name instead of silently falling back to `balanced`; `ScoringConfig` rejects negative weights and `gc_min > gc_max` (library-direct callers only — no product-path behavior change)
- Guarded `calculate_mfe()` against unbounded ViennaRNA `RNA.fold()` runtime on long sequences (algorithmic-complexity DoS); sequences over 1000nt now skip MFE folding instead of risking an indefinite worker hang (correction: this also shipped in v3.2.7 but was undocumented here at release time)

### Added
- Benchmark runs now record whether ViennaRNA bindings were active, going forward

## v3.2.6 — 2026-06-27

### Changed
- Aligned the public-safe wet-lab feedback intake surface with the actual GitHub issue template, removing fields that were never collected
- Disabled blank GitHub issues to reduce the chance of sensitive wet-lab/construct details landing in an untemplated public issue

## v3.2.5 — 2026-06-25

### Fixed
- The published Docker image's `POST /api/optimize` crashed with an `AttributeError` due to a routing bug in the local dev/Docker server (hosted web app and PyPI CLI were unaffected)
- `CITATION.cff`'s `doi` field corrected to the v3.2.4 exact-release DOI

## v3.2.4 — 2026-06-24

### Added
- Runtime validation registry exposing all 17 sequence checks (9 advisory scanners + restriction-site + MoClo-overhang) via the web "Sequence Checks" panel and `/api/optimize` response
- Explicit consent modal before sending sequences to AlphaFold DB / ESM Atlas structure-prediction links

### Fixed
- Web UI "AA Preserved" badge showed assembly/restriction-site review status instead of actual amino acid identity
- `/api/optimize/compare` and `/api/optimize/batch` now reject `host`/`host_profile` fields with HTTP 400 instead of silently dropping them
- `high_cai` profile and `feasibility_best` objective are now consistently disclosed/rejected for non-default hosts across CLI, REST, and web UI

## v3.2.3 — 2026-06-19

### Fixed
- Release provenance hashing now computed from the committed git blob instead of local working-tree bytes, fixing CRLF/LF drift on Windows that could silently produce incorrect SHA-256 values in reproducibility manifests
- Public-surface DOI references switched from version-pinned Zenodo DOIs to the concept DOI, which always resolves to the latest release

### Added
- A public-surface audit and a CHANGELOG consistency check now run on every push/PR in CI

## v3.2.2 — 2026-06-18

### Fixed
- **`multi_constraint_pass` definition corrected (scoring_contract v1.1)** — now requires GC-in-target-range in addition to biological/assembly pass; corrected ablation values (L3=3.5%, L4=5.6%) supersede the previous inflated figures. All benchmark artifacts regenerated (seed=320, N=49,257); Zenodo `benchmark_results.csv` v2 (DOI: [10.5281/zenodo.20676276](https://doi.org/10.5281/zenodo.20676276)) supersedes v1
- Benchmark: source-profile codon-table injection now flows into both the design and scoring paths, fixing a prior gap where design always used the default table regardless of an injected profile

### Added
- Data: three genome-annotated *N. benthamiana* codon-usage profiles (SGN QLD183 v103 CDS-derived; SGN NbeV1.1 all-CDS-derived; SGN NbeV1.1 high-confidence-CDS-derived), alongside the existing packaged reference profile

### Changed
- Web: Host System cards now render dynamically from the API's `supported_hosts`/`host_metadata` instead of being hardcoded

### Docs
- Aligned public README/docs/web/citation/packaging/roadmap/benchmark wording with the in-silico CDS design claim boundary
- Clarified that public GitHub Issues must not contain raw sequences, confidential construct details, or other sensitive data

---

## v3.2.1 — 2026-06-16

### Added
- **Protein risk annotation** — transmembrane helix (Kyte-Doolittle, window=19) and signal peptide heuristics flag HIGH / MEDIUM / LOW / UNKNOWN risk sequences before wet-lab submission

### Fixed
- CAI provenance annotation in benchmark output
- Type IIS restriction site warning status
- Manifest SHA-256 reproducibility drift on Windows
- wet-lab result GitHub template: protein_class options + validation consent checkbox

---

## v3.2.0 — 2026-06-11

### Benchmark
- **Formal benchmark** — N. benthamiana SGN CDS dataset (N=49,257 sequences), seed=320 for deterministic reproduction. All metrics are in-silico; no wet-lab validation is claimed.
- **Seed injection** — benchmark runner now accepts `--seed` for fully deterministic reruns; `most_frequent_codon` tie-breaking deduplication included

### Data / Provenance
- **Codon table provenance disclosure** — `nbenthamiana_codons.json` origin documented as a legacy Kazusa + SGN v1.0.1-era reference; `codon_table_manifest.json` added with sha256 pin, `build_path_status: incomplete`, and known limitations recorded

### API / Schema
- **Design Package schema v1.0.0** — formal IUPAC/FASTA I/O contracts and MFE null invariant established
- **Registry constants export** — `DEFAULT_CAI_TARGET`, `DEFAULT_GC_LOW`, `DEFAULT_GC_HIGH` importable as public production constants

### Fixed
- **Input validator** — IUPAC ambiguous DNA/AA sequence misclassification corrected

### Documentation
- **Benchmark Trust Pack** — registry/spec integrity tests, metric correctness tests, and deterministic baseline runner added
- **Claim wording alignment** — public-facing API and CLI output wording unified; no expression-level or yield improvement claims

---

## v3.1.9 — 2026-06-04

### Documentation
- Internal housekeeping: project tracking references updated. No engine changes.

---

## v3.1.8 — 2026-06-01

### Breaking
- **`gc_target` default** — now targets host GC midpoint (~60%) instead of legacy 42.5%

### Changed
- GC scoring replaced with band function; `assembly_friendly` weights updated
- `feasibility.py` defaults realigned to 55–65% output distribution

### Fixed
- CAI reference unified to golden-set weights throughout
- Duplicate root `data/` files removed

## v3.1.7 — 2026-05-31

### Added
- **Host selector UI** — N. benthamiana / BY-2 Experimental toggle; BY-2 disables Feasibility Best
- **E2E smoke tests** — 5 Playwright tests auto-run after each deployment
- **Eijex MCP** — MCP-compatible access via mcp.eijex.com

## v3.1.6 — 2026-05-30

### Added
- **BY-2 host (experimental)** — `--host by2`, N. tabacum Kazusa codon table
- **Profile comparison** — `--compare-profiles` CLI + `POST /api/optimize/compare`
- **Batch API** — `POST /api/optimize/batch` (최대 20개 서열)
- **Structure links** — AlphaFold DB + ESM Atlas 링크 결과 화면에 추가
- **Experimental scoring hooks** — non-public optional scoring hooks; public profiles unchanged
- **GFP Tutorial** — `docs/tutorials/gfp-nbenthamiana.md`

### Fixed
- API GC 기본값 CLI와 일치 (40-55% → 55-65%)

## v3.1.5 — 2026-05-28

### Fixed
- **Data file packaging** — PyPI wheel에 JSON 데이터 파일 미포함 버그 수정 (`pip install` 후 경로 오류)
- **Path resolution** — `parents[N]` 하드코딩 → `get_data_path()` 통일 (pip/Docker/Vercel 모두 정상 작동)
- **CI matrix** — Python 3.13 + Windows 테스트 추가

## v3.1.4 — 2026-05-27

### Added
- **CITATION.cff** — GitHub "Cite this repository" button; version-bumped automatically
- **SECURITY.md** — vulnerability reporting policy
- **ROADMAP.md** — public development direction, validation scope, and planned host/profile work
- **scripts/release.py** — automates version string updates

### Changed
- **Public engine naming** — `factorforge.engines.profile` / CLI `--engine profile` formalized
- **Benchmark corrected** — CAI 0.76, GC% 59.77%, target 55–65% (N=3,876 SGN CDS, balanced profile)

### Fixed
- **Stale test version** — `test_design_package.py` product_version 3.1.1 → 3.1.3

## v3.1.3 — 2026-05-26

### Fixed
- **Disabled profile cards** — tooltips now accessible; "Pending wet-lab validation" notice added to 5' Ramp and Viral Delivery
- **Viral Delivery tooltip** — corrected citation reference; updated to Peccoud et al. 2024 (PMC11718241)
- **Analytics notice** — clarified: "Submitted sequences are not logged or stored"

## v3.1.2 — 2026-05-26

### Fixed
- **viral_delivery scoring** — corrected citation reference; `w_mfe` 0.40→0.30 per PMC11718241 (Peccoud et al. 2024)
- **5' Ramp deoptimization** — N-terminal deoptimization bottom 50%→25% (mild, per PMC11718241)

## v3.1.1 — 2026-05-24

### Added
- **Wet-lab feedback link** — Submit result button opens a public-safe feedback path
- **JSON Copy button** — one-click copy of full optimization JSON output

### Changed
- **Design Objective order** — reordered to match recommended wet-lab testing sequence
- **Validation fields** — Issue template updated for structured public-safe feedback

### Fixed
- **Vercel deployment** — resolved /api/optimize 404 caused by incorrect Root Directory setting

## v3.1.0 — 2026-05-24

### Added
- **Custom restriction site removal** — synonymous substitution of user-specified restriction sites
- **Rare codon run detection** — detects consecutive rare codons for ribosome stalling risk
- **Dinucleotide reduction** — CAI-budgeted CpG/TpA reduction (aggressive / balanced / cai_preserving modes)

### Fixed
- **Pipeline metric accuracy** — CAI/GC/score re-measured from final sequence after dinucleotide fix
- **CAI guard weight consistency** — guard now uses Sharp & Li 1987 golden set reference weights

## v3.0.0 — 2026-05-23

First official release. DP Feasibility Engine, CLI, Validator, open-source under AGPL-3.0.
