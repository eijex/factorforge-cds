# Changelog

Full changelog: [CHANGELOG.md on GitHub](https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md)

FactorForge v3.0+ is the current public release line.

## v3.2.1 — 2026-06-16

### Added
- **Protein risk annotation** — transmembrane helix (Kyte-Doolittle, window=19) and signal peptide heuristics flag HIGH / MEDIUM / LOW / UNKNOWN risk sequences before wet-lab submission

### Fixed
- CAI provenance annotation in benchmark output
- Type IIS restriction site warning status
- Manifest SHA-256 reproducibility drift on Windows
- wet-lab result GitHub template: protein_class options + validation consent checkbox

### Changed
- Google Form wet-lab submission channel added

---

## v3.2.0 — 2026-06-11

### Benchmark
- **Formal benchmark** — N. benthamiana SGN CDS dataset (N=49,257 sequences), seed=320 for deterministic reproduction. All metrics are in-silico; no wet-lab validation is claimed.
- **Seed injection** — benchmark runner now accepts `--seed` for fully deterministic reruns; `most_frequent_codon` tie-breaking deduplication included (099)

### Data / Provenance
- **Codon table provenance disclosure** — `nbenthamiana_codons.json` origin documented as a legacy Kazusa + SGN v1.0.1-era reference; `codon_table_manifest.json` added with sha256 pin, `build_path_status: incomplete`, and known limitations recorded (097)

### API / Schema
- **Design Package schema v1.0.0** — formal IUPAC/FASTA I/O contracts and MFE null invariant established (090)
- **Registry constants export** — `DEFAULT_CAI_TARGET`, `DEFAULT_GC_LOW`, `DEFAULT_GC_HIGH` importable as public production constants (091)

### Fixed
- **Input validator** — IUPAC ambiguous DNA/AA sequence misclassification corrected (098)

### Documentation
- **Benchmark Trust Pack** — registry/spec integrity tests, metric correctness tests, and deterministic baseline runner added (084)
- **Claim wording alignment** — public-facing API and CLI output wording unified; no expression-level or yield improvement claims (092)

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
