# Roadmap

FactorForge development direction and planned work.

> Current release: **v3.1.5** — [Changelog](CHANGELOG.md) · [Releases](https://github.com/eijex/factorforge-cds/releases)
> GitHub Milestones: [github.com/eijex/factorforge-cds/milestones](https://github.com/eijex/factorforge-cds/milestones)

---

## Current Focus (v3.1.x)

Wet-lab validation and stability. The core engine is feature-complete; this track collects experimental results and applies targeted fixes based on real expression data.

- Wet-lab feedback integration (Google Form → GitHub Issues)
- Bug fixes and metric corrections based on community reports
- Documentation improvements

---

## v3.2 — Accessibility + Foundation

| Feature | Description |
|---------|-------------|
| **Bioconda package** | Conda install support for bioinformatics users |
| **Batch optimization API** | Multi-sequence input for pipeline use cases |
| **Profile comparison mode** | Side-by-side results across multiple profiles for the same input |
| **Additional plant hosts** | BY-2 cells + further host additions based on expression system research |
| **Tutorial / worked example** | At least one end-to-end worked example in documentation |
| **CODE_OF_CONDUCT.md** | Community guidelines |

---

## v3.3 — Algorithm Depth

| Feature | Description |
|---------|-------------|
| **tAI support** | tRNA Adaptation Index as alternative/complement to CAI |
| **5' UTR mRNA folding** | Active optimization of secondary structure around start codon |
| **Codon pair bias** | Detection and adjustment of unfavorable adjacent codon pairs |
| **Kozak context** | Optimization of nucleotides flanking the ATG (-3, +4 positions) |
| **Extended DP engine** | Expose DP feasibility parameters via CLI and API |

---

## v3.4 — Protein Type Expansion

| Feature | Description |
|---------|-------------|
| **Membrane protein profile** | TMD detection; reduced translation speed at transmembrane domains |
| **Secreted protein profile** | Signal peptide recognition; separate optimization of SP region |
| **Chloroplast-targeted profile** | Transit peptide-aware optimization |
| **PTGS silencing risk** | Flag sequences with similarity to plant endogenous genes |

---

## v3.5 — Validation-Driven

| Feature | Description |
|---------|-------------|
| **5' Ramp activation** | Enable once wet-lab data confirms N-terminal ramp benefit |
| **Viral Delivery activation** | Enable once wet-lab data confirms scoring validity |
| **Tissue-specific optimization** | Leaf vs seed vs root codon preference tuning |
| **Publication draft** | Software manuscript draft |

---

## v3.6 — Host Expansion

| Feature | Description |
|---------|-------------|
| **Insect expression (Sf9, Tni)** | Baculovirus/BEVS codon optimization |
| **Additional plant hosts** | Arabidopsis, tomato, Lemna, Wolffia globosa |
| **Publication finalized** | Manuscript finalized |

---

## v3.7 — ML Engine *(data-conditional)*

A ML-based CDS design engine added as `--engine ml` once sufficient wet-lab validation data is available. The v3-alpha ML prototype (archived under `archive/v3-ml-prototype/`) informs the architecture.

- ML engine exposed as `--engine ml` alongside the existing rule-based engine
- Training pipeline based on accumulated wet-lab validation data
- Rule-based vs ML benchmark comparison
- No committed timeline — depends on wet-lab data accumulation

---

## v3.8 — Release Readiness

- Full API docstring coverage
- Reproducible benchmark scripts
- Extended tutorial
- State of the field comparison (JCat, OPTIMIZER, Codon Harmonizer, IDT)

---

## v3.9 — Stable Release

- All profiles wet-lab validated and active
- Full documentation audit
- Final polish

---

## Wet-lab Validation *(ongoing — not version-bound)*

Wet-lab validation is an ongoing, open-ended process independent of software versioning.
Results are tracked as GitHub Issues under the [Wet-lab Validation](https://github.com/eijex/factorforge-cds/milestone/9) milestone.

- Submit results via [Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform?usp=header) or [GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml)
- See [VALIDATION.md](VALIDATION.md) for submission format and current status

---

## Out of Scope

- Animal or clinical expression systems
- Protein structure prediction (use AlphaFold / ESMFold directly)
- Codon optimization for microbial hosts (E. coli, yeast)
