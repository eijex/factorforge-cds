# Validation

## Status

FactorForge is an **in-silico only** tool. It optimizes codon usage, GC content, and sequence constraints — but it cannot guarantee protein expression, yield, solubility, or folding in your specific system.

!!! warning
    Wet-lab validation is always required before use in production workflows.

## What FactorForge Checks

| Check | Covered | Notes |
|-------|---------|-------|
| Internal stop codons | ✅ | Hard fail |
| Amino acid identity | ✅ | Hard fail |
| GC% range | ✅ | Host-dependent; default *N. benthamiana*: 55–65% |
| Forbidden restriction sites | ✅ | BsaI, BsmBI, BpII (Golden Gate) |
| Invalid codons | ✅ | Hard fail |
| Rare codon runs | ✅ | Ribosome stalling risk detection |
| CpG/TpA dinucleotides | ✅ | CAI-budgeted reduction |
| Protein folding | ❌ | Requires wet-lab or AlphaFold |
| Actual expression level | ❌ | Requires wet-lab |
| Yield / solubility | ❌ | Requires wet-lab |
| Signal peptide behavior | ❌ | Requires wet-lab |
| Membrane protein topology | ❌ | Requires wet-lab |

## Contribute Wet-Lab Results

If you tested a FactorForge-optimized sequence in the lab, please share your result — positive or negative. Failed experiments are just as valuable.

**→ [Submit via Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform?usp=header)** (recommended)

Or open a [GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) with the label `wet-lab-result`.

Please include:

- Host organism and protein
- FactorForge version and optimization profile
- Expression system and assay
- Result (qualitative or quantitative)
- Promoter, subcellular targeting, harvest timepoint (optional)
- Whether a native/unoptimized sequence was used as a control

Contributors will be credited in [VALIDATION.md](https://github.com/eijex/factorforge-cds/blob/main/VALIDATION.md) and future releases.
