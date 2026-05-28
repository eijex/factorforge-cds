# Changelog

Full changelog: [CHANGELOG.md on GitHub](https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md)

FactorForge v3.0+ is the public release of an internal deterministic rule-based engine (v2) that succeeded an earlier thesis-derived prototype (v1) and a later ML prototype (v3-alpha, archived). See [README](https://github.com/eijex/factorforge-cds#development-history) for the full version lineage.

## v3.1.4 — 2026-05-27

### Added
- **CITATION.cff** — GitHub "Cite this repository" button; version-bumped automatically
- **SECURITY.md** — vulnerability reporting policy
- **ROADMAP.md** — public development direction (v3.2–v3.9, ML engine v3.7+, wet-lab scope)
- **bump_version.py** — automates version string updates across 14 files
- **Development history** — v1→v2→v3-alpha→v3.x→v3.7+(ML) lineage documented across README, CHANGELOG, archive READMEs

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
- **Wet-lab feedback modal** — Submit result button opens embedded Google Form with version and profile pre-filled
- **JSON Copy button** — one-click copy of full optimization JSON output

### Changed
- **Design Objective order** — reordered to match recommended wet-lab testing sequence
- **Validation fields** — Issue template updated with promoter, subcellular targeting, harvest timepoint, native control

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
