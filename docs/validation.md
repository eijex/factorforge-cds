# Validation

## Status

FactorForge is an **in-silico only** CDS design tool. It optimizes codon usage, GC content, and sequence constraints, but it cannot guarantee protein expression, yield, solubility, folding, secretion, or biological activity in any specific system.

!!! warning
    Wet-lab validation is always required before use in production workflows.

Public validation entries are manually reviewed, non-confidential wet-lab feedback summaries. They are not controlled validation studies, regulatory claims, or guarantees of FactorForge performance.

## What FactorForge Checks

| Check | Covered | Notes |
|-------|---------|-------|
| Internal stop codons | Yes | Hard fail |
| Amino acid identity | Yes | Hard fail |
| GC% range | Yes | Host-dependent; default *N. benthamiana*: 55-65% |
| Forbidden restriction sites | Yes | BsaI, BsmBI, BpiI (Golden Gate) |
| Invalid codons | Yes | Hard fail |
| Rare codon runs | Yes | Ribosome stalling risk detection |
| CpG/TpA dinucleotides | Yes | CAI-budgeted reduction |
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
- Expression system / assay type.
- Harvest timepoint.
- Replicate category.
- Comparison result.
- Expression result.
- Approximate yield range.
- Non-confidential notes only.

Submission links:

- [Submit via Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform) — structured wet-lab result submission
- [Open a public GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) — public-safe coarse summaries only
- Email `eijex.lab@gmail.com` — private or sensitive summaries
- See [VALIDATION.md](https://github.com/eijex/factorforge-cds/blob/main/VALIDATION.md) for the current public validation log.

All public entries require manual review before publication.
