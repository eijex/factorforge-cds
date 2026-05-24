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
| GC% range | ✅ | Target 40–60% for *N. benthamiana* |
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

- Open a [GitHub Issue](https://github.com/eijex/factorforge-cds/issues) with the label `wet-lab-result`
- Or email: eijex.lab@gmail.com

**Include:**

- Host organism
- Protein name / class (e.g. antibody, enzyme, membrane protein)
- FactorForge version and run parameters
- Expression system and assay used
- Result (qualitative or quantitative)

Contributors will be credited in [VALIDATION.md](https://github.com/eijex/factorforge-cds/blob/main/VALIDATION.md) and future releases.
