# Evidence Intake Policy

This document defines what counts as acceptable evidence for FactorForge benchmark and validation claims, and how new evidence is incorporated.

## Acceptable evidence sources

- Public, peer-reviewed literature with documented methodology
- Public sequence databases and reference datasets with stable accessions
- Reproducible computational benchmarks committed to this repository (scripts, fixtures, registry-backed configuration)
- Wet-lab results that are documented with methodology, raw data availability, and attribution

## Intake process

1. Evidence is reviewed for source credibility, reproducibility, and relevance to the specific claim it would support.
2. Computational evidence must be reproducible from committed scripts and either committed fixtures or documented public data sources.
3. Evidence that supports a public claim is linked from the relevant documentation (for example `benchmarks/README.md` or `docs/benchmark.md`).
4. Evidence that does not meet these standards is not used to support public claims, even if it informs internal development decisions.

## Non-acceptable evidence for public claims

- Private, non-reproducible internal results
- Anecdotal or informal observations without documented methodology
- Third-party claims without verifiable sourcing

## Relationship to the Public Claim Policy

This policy determines what evidence is admissible. [`benchmarks/PUBLIC_CLAIM_POLICY.md`](../../benchmarks/PUBLIC_CLAIM_POLICY.md) determines how admissible evidence may be represented publicly. Both apply together when evaluating any new benchmark or performance claim.
