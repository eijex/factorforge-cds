# Validation

## Status

FactorForge is an **in-silico only** CDS design tool. It optimizes codon usage, GC content, and sequence constraints, but it cannot guarantee protein expression, yield, solubility, folding, secretion, or biological activity in any specific system.

!!! warning
    Wet-lab validation is always required before use in production workflows.

Public validation entries are manually reviewed, non-confidential wet-lab feedback summaries. They are not controlled validation studies, regulatory claims, or guarantees of FactorForge performance.

## What FactorForge Checks

Checks are grouped by domain (see `rule-engine-roadmap.md` for full per-rule
status).

**Sequence integrity** (hard fail in the production design pipeline):

| Check | Covered | Notes |
|-------|---------|-------|
| Internal stop codons | Yes | Hard fail |
| Amino acid identity | Yes | Hard fail |
| Invalid / partial codons | Yes | Hard fail |
| Reading-frame / length consistency | Yes | Hard fail |

**Configured constraints**:

| Check | Covered | Notes |
|-------|---------|-------|
| GC% range | Yes | Host-dependent; default *N. benthamiana*: 55-65% (legacy engine-output-calibrated band, `v1`; see Codon Reference Contract below). *N. tabacum* (BY-2): 55-65% (unchanged). This range is a configured target metric, not a hard gate outside the benchmark scoring contract. |

### Codon Reference Contract (v1 / v2)

`codon_reference_contract_version` tracks which **codon usage table and GC
reference band** a result was generated against — distinct from
`scoring_contract_version` (currently `v1.1`), which defines the unrelated
`multi_constraint_pass` pass/fail formula (`biological_pass AND assembly_pass
AND gc_in_target_range`). A v1→v2 codon-reference change does not change what
"pass" means; it changes the codon table and GC band a result is scored
against.

| Contract version | Codon reference asset | GC reference band (*N. benthamiana*) | Status |
|---|---|---|---|
| `v1` | `nbenthamiana_legacy_kazusa_sgn_v101` (legacy, circular-derived) | 55-65% | historical comparator retained for v3.2.x continuity and sensitivity interpretation; not the current default |
| `v2` | `nbenthamiana_nbev11_hc_v2` (NbeV1.1 LAB-strain, native genome-composition anchor; `reference-policy audit`) | 40-47% | current production software default; public overrides remain disabled; not experimental validation or comparative biological-performance evidence |

The active default is recorded in `data/reference/active_codon_reference.json`
and synchronized with `current_parameter_registry.yaml`'s `codon_reference.active`
block (see `tests/test_registry_production_sync.py` for the sync guard). Public
CLI and REST API optimize endpoints do **not** expose a codon-reference override;
they always use the active default and reject request fields such as
`codon_reference`, `codon_table_id`, or `codon_table_path`. Python/benchmark
research paths can still pass an explicit `codon_table_path` for controlled
comparator runs, but those packaged assets are not public production options
unless a separate activation gate enables them.

**Advisory sequence-risk scans** (all 9 run by default; findings are reported, never gating):

| Check | Covered | Notes |
|-------|---------|-------|
| PolyA-like motifs | Yes | Heuristic |
| AU-rich elements (ARE) | Yes | Heuristic |
| AT-rich runs | Yes | Minimum run length 6 nt default |
| Homopolymers | Yes | Synthesis-risk threshold |
| Tandem repeats | Yes | >=15 nt, recombination risk |
| Local GC extremes | Yes | 50-nt window, 25-75% synthesis guard |
| Splice-like motifs | Yes | Heuristic; false positive risk noted |
| CpG/TpA dinucleotides | Yes | CAI-budgeted reduction |
| Rare codon runs | Yes | Ribosome stalling risk detection |

**Assembly review**:

| Check | Covered | Notes |
|-------|---------|-------|
| Forbidden restriction sites | Yes | BsaI, BsmBI, BpiI (Golden Gate); halts the production design pipeline when unresolvable |
| MoClo overhang validity / collision | Yes (opt-in) | Not run by default; reports warnings only, never halts the pipeline |

**Outside current scope** (requires wet-lab validation):

| Check | Covered | Notes |
|-------|---------|-------|
| Protein folding | No | Requires wet-lab validation |
| Actual expression level | No | Requires wet-lab validation |
| Yield / solubility | No | Requires wet-lab validation |
| Signal peptide behavior | No | Requires wet-lab validation |
| Membrane protein topology | No | Requires wet-lab validation |

## Contribute Wet-Lab Results

Public-safe validation summaries are welcome. Positive, negative, equivalent, reduced, failed-expression, and inconclusive outcomes are all useful.

Public GitHub Issues are public. Use them only for public-safe summaries. Do not submit private or sensitive wet-lab feedback through public GitHub Issues.

Default public fields should use protein class rather than protein name. Public credit requires explicit approval and maintainer review.

Do not submit:

- Raw DNA, RNA, or protein sequences.
- Proprietary or unpublished construct sequences.
- Internal batch IDs.
- Confidential partner/customer information.
- Patient or clinical subject data.
- Private contact information.
- Exact confidential process parameters.
- Confidential assay protocols.
- Unpublished construct-identifying labels.

Public-safe summaries may include:

- FactorForge version.
- Host organism.
- Optimization profile.
- Protein class.
- Coarse assay category.
- Replicate category.
- Comparison result.
- Expression result.
- Non-confidential notes only.

Submission channels:

- [Share Wet-lab Results (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) — public-safe coarse summaries only
- Email `eijex.lab@gmail.com` — private or sensitive summaries
- See [VALIDATION.md](https://github.com/eijex/factorforge-cds/blob/main/VALIDATION.md) for the current public validation log.

All public entries require manual review before publication.
