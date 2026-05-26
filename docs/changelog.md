# Changelog

Full changelog: [CHANGELOG.md on GitHub](https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md)

## v3.1.3 — 2026-05-26

### Fixed
- **Disabled profile cards** — tooltips now accessible; "Pending wet-lab validation" notice added to 5' Ramp and Viral Delivery
- **Viral Delivery tooltip** — corrected citation reference; updated to Peccoud et al. 2024 (PMC11718241)
- **Analytics notice** — clarified: "Submitted sequences are not logged or stored"

## v3.1.2 — 2026-05-26

### Fixed
- **viral_delivery scoring** — corrected citation reference; `w_mfe` 0.40→0.30 per PMC11718241 (Peccoud et al. 2024)
- **5' Ramp deoptimization** — N-terminal deoptimization bottom 50%→25% (mild, per PMC11718241)

### Changed
- **Changelog label** — "ML Research Track" renamed to "Research Track"

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
