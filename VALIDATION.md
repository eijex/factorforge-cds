# Validation Log

FactorForge is an in-silico CDS design tool. It optimizes codon usage, GC content, and sequence constraints — but it cannot guarantee protein expression, yield, solubility, or folding in your specific system. Wet-lab validation is always required.

This file tracks community wet-lab results to build a shared evidence base over time.

**Try it first**: [factorforge-cds.vercel.app](https://factorforge-cds.vercel.app)

---

## What FactorForge checks (and what it doesn't)

| Check | Covered | Notes |
|-------|---------|-------|
| Internal stop codons | ✅ | Hard fail |
| Amino acid identity | ✅ | Hard fail |
| GC% range | ✅ | Target 40–60% for *N. benthamiana* |
| Forbidden restriction sites | ✅ | BsaI, BsmBI, BpII (Golden Gate) |
| Invalid codons | ✅ | Hard fail |
| Protein folding | ❌ | Requires wet-lab or AlphaFold |
| Actual expression level | ❌ | Requires wet-lab |
| Yield / solubility | ❌ | Requires wet-lab |
| Signal peptide behavior | ❌ | Requires wet-lab |
| Membrane protein topology | ❌ | Requires wet-lab |

Passing all checks means the sequence is **biologically plausible and assembly-ready** — not that it will express well.

---

## How to Contribute

If you tested a FactorForge-optimized sequence in the lab, please share your result — positive or negative. Failed experiments are just as valuable.

- Open a [GitHub Issue](https://github.com/eijex/factorforge-cds/issues) with the label `wet-lab-result`
- Or email: eijex.lab@gmail.com

Include:
- Host organism
- Protein name / class (e.g. antibody, enzyme, membrane protein)
- FactorForge run parameters or optimized sequence
- Expression system and assay used
- Result (qualitative or quantitative)

Contributors will be credited below.

---

## Results

*No results yet. Be the first to validate FactorForge in the lab.*

---

## Credits

| Contributor | Institution | Protein | Host | Result | Date |
|-------------|-------------|---------|------|--------|------|
| — | — | — | — | — | — |
