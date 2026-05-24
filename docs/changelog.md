# Changelog

Full changelog: [CHANGELOG.md on GitHub](https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md)

## v3.1.0 — 2026-05-24

### Added
- **Custom restriction site removal** — synonymous substitution of user-specified restriction sites
- **Rare codon run detection** — detects consecutive rare codons for ribosome stalling risk
- **Dinucleotide reduction** — CAI-budgeted CpG/TpA reduction (aggressive / balanced / cai_preserving modes)

### Fixed
- **Pipeline metric accuracy** — CAI/GC/score re-measured from final sequence after dinucleotide fix
- **CAI guard weight consistency** — guard now uses Sharp & Li 1987 golden set reference weights

## v3.0.0 — 2026-05-23

First official release. DP Feasibility Engine, CLI, Validator, open-source under Apache 2.0.
