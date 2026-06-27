# Validation Log

FactorForge is an in-silico CDS design tool. It optimizes codon usage, GC content, and sequence constraints, but it cannot guarantee protein expression, yield, solubility, folding, secretion, or biological activity in any specific system. Wet-lab validation is always required.

This file tracks manually reviewed, non-confidential wet-lab feedback summaries. Public entries are user-submitted feedback summaries, not controlled validation studies, regulatory claims, or guarantees of FactorForge performance.

**Try it first**: [factorforge.eijex.com](https://factorforge.eijex.com)

---

## What FactorForge Checks

| Check | Covered | Notes |
|-------|---------|-------|
| Internal stop codons | Yes | Hard fail |
| Amino acid identity | Yes | Hard fail |
| GC% range | Yes | Host-dependent; default *N. benthamiana*: 40-47% (native genome-composition reference band, v3.3.0+, not a wet-lab validated expression optimum). *N. tabacum* (BY-2): 55-65% (unchanged). |
| Forbidden restriction sites | Yes | BsaI, BsmBI, BpiI (Golden Gate) |
| Invalid codons | Yes | Hard fail |
| Protein folding | No | Requires wet-lab validation |
| Actual expression level | No | Requires wet-lab validation |
| Yield / solubility | No | Requires wet-lab validation |
| Signal peptide behavior | No | Requires wet-lab validation |
| Membrane protein topology | No | Requires wet-lab validation |

Passing all checks means the sequence is biologically plausible and assembly-ready. It does not mean the sequence will express well.

---

## Public-Safe Contributions

Public-safe validation summaries are welcome. Positive, negative, equivalent, reduced, failed-expression, and inconclusive outcomes are all useful when described without confidential details.

Public GitHub Issues are public. Use them only for public-safe summaries. Do not submit private or sensitive wet-lab feedback through public GitHub Issues.

If using email, submit only information you have permission to share and do not include raw sequences or private identifiers unless there is an explicit private review agreement.

## Do Not Submit

- Raw DNA, RNA, or protein sequences.
- Proprietary or unpublished construct sequences.
- Internal batch IDs.
- Confidential partner/customer information.
- Patient or clinical subject data.
- Private contact information.
- Private customer names.
- Exact confidential process parameters.
- Confidential assay protocols.
- Unpublished construct-identifying labels.

## Public-Safe Fields

- FactorForge version.
- Host organism.
- Optimization profile.
- Protein class, not confidential protein identity.
- Expression system / assay type.
- Harvest timepoint.
- Replicate category.
- Comparison result.
- Expression result.
- Approximate yield range.
- Non-confidential notes only.

Contributor identity is not listed in public summaries by default. Public credit requires explicit approval and additional maintainer review.

Submission channels:

- [Share Wet-lab Results (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) — public-safe coarse summaries only
- Email `eijex.lab@gmail.com` for private or sensitive summaries

---

## Manual Review Checklist

- Optional public-summary consent is present.
- No raw DNA/RNA/protein sequence is present.
- No internal batch ID is present.
- No confidential construct identifier is present.
- No customer or partner name is present without permission.
- No patient or clinical subject data is present.
- No submitter name, respondent email, phone number, or contact information is present.
- No exact confidential process parameter is present.
- Protein name is removed or coarse-grained unless explicitly approved after maintainer review.
- Yield is expressed as a coarse range.
- Text does not imply controlled validation, regulatory approval, or guaranteed performance.
- Disclaimer language is preserved.

---

## Results

No manually reviewed public-safe validation summaries have been published yet.

---

## Public Summary Disclaimer

Public validation summaries are user-submitted wet-lab feedback summaries. They are not controlled validation studies, regulatory claims, or guarantees of FactorForge performance.
