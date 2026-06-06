# Validation Log

FactorForge is an in-silico CDS design tool. It optimizes codon usage, GC content, and sequence constraints — but it cannot guarantee protein expression, yield, solubility, or folding in your specific system. Wet-lab validation is always required.

This file tracks community wet-lab results to build a shared evidence base over time.

**Data use notice:** Submitted results may be used in aggregate to improve FactorForge. No personally identifying information is collected. Sequences are not stored or logged.

**Try it first**: [factorforge.eijex.com](https://factorforge.eijex.com)

---

## What FactorForge checks (and what it doesn't)

| Check | Covered | Notes |
|-------|---------|-------|
| Internal stop codons | ✅ | Hard fail |
| Amino acid identity | ✅ | Hard fail |
| GC% range | ✅ | Host-dependent; default *N. benthamiana*: 55–65% |
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

- **[Submit via Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform?usp=header)** (recommended)
- Or open a [GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) with the label `wet-lab`
- Or email: eijex.lab@gmail.com

Include:
- Host organism
- Protein name / class (e.g. antibody, enzyme, membrane protein)
- FactorForge version and optimization profile used
- Expression system and assay used
- Result (qualitative or quantitative)
- Promoter used (optional)
- Subcellular targeting (optional)
- Harvest timepoint (optional)
- Whether a native/unoptimized sequence was used as a control

Contributors will be credited below.

---

## Results

*No results yet. Be the first to validate FactorForge in the lab.*

---

## Credits

| Contributor | Institution | Protein | Host | Result | Date |
|-------------|-------------|---------|------|--------|------|
| — | — | — | — | — | — |
