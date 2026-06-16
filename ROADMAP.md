# Roadmap

FactorForge development direction and planned work.

> Current release: **v3.2.1** — [Changelog](CHANGELOG.md) · [Releases](https://github.com/eijex/factorforge-cds/releases)
> GitHub Milestones: [github.com/eijex/factorforge-cds/milestones](https://github.com/eijex/factorforge-cds/milestones)

---

## Current Focus (v3.2.x)

Public-claim alignment, wet-lab feedback intake, and stability. The core engine is feature-complete for the current public scope; this track keeps public documentation, package metadata, web UI, Docker/GHCR references, and citation archives aligned with the same in-silico CDS design boundary.

- Public surface wording audit across README, docs, web UI, package metadata, Docker/GHCR references, and citation archives
- Wet-lab feedback review and public-safe validation summaries
- Raw-sequence logging and benchmark-summary leakage checks
- Codon-usage table provenance notes and checksums for public datasets
- Bug fixes and metric corrections based on community reports
- Documentation improvements

---

## v3.2 — Accessibility + Foundation

| Feature | Description |
|---------|-------------|
| **Bioconda package** | Conda install support for bioinformatics users |
| **Batch optimization API** | Multi-sequence input for pipeline use cases |
| **Profile comparison mode** | Side-by-side results across multiple profiles for the same input |
| **Additional plant hosts** | BY-2 cells + further host additions based on host-context and codon-source research |
| **Tutorial / worked example** | At least one end-to-end worked example in documentation |
| **CODE_OF_CONDUCT.md** | Community guidelines |
| **Pre-submission hardening checklist** | Version-aligned public surfaces, benchmark reproducibility, security wording, and codon-table provenance |

---

## v3.3 — Algorithm Depth

| Feature | Description |
|---------|-------------|
| **tAI support** | tRNA Adaptation Index as alternative/complement to CAI |
| **5' UTR mRNA folding** | Optional MFE-aware review around the start codon when dependencies are available |
| **Codon pair bias** | Detection and adjustment of unfavorable adjacent codon pairs |
| **Kozak context** | Optimization of nucleotides flanking the ATG (-3, +4 positions) |
| **Extended DP engine** | Expose DP feasibility parameters via CLI and API |

---

## v3.4 — Protein Type Expansion

| Feature | Description |
|---------|-------------|
| **Membrane protein profile** | TMD detection; local codon-usage review around transmembrane domains |
| **Secreted protein profile** | Signal peptide recognition; separate optimization of SP region |
| **Chloroplast-targeted profile** | Transit peptide-aware optimization |
| **PTGS silencing risk** | Flag sequences with similarity to plant endogenous genes |

---

## v3.5 — Validation-Driven

| Feature | Description |
|---------|-------------|
| **5' Ramp activation** | Enable only if reviewed wet-lab feedback supports N-terminal ramp benefit |
| **Viral Delivery activation** | Enable only if reviewed wet-lab feedback supports scoring validity |
| **Tissue-specific design profiles** | Leaf vs seed vs root codon-source review and profile gating |
| **Documentation audit** | Comprehensive docs and capability review |

---

## v3.6 — Host Expansion

| Feature | Description |
|---------|-------------|
| **Insect host profiles (Sf9, Tni)** | Baculovirus/BEVS host-context research and codon-source review |
| **Additional plant hosts** | Arabidopsis, tomato, Lemna, Wolffia globosa |
| **Documentation finalized** | Final docs pass |

---

## v3.7 — Release Readiness

- Full API docstring coverage
- Reproducible benchmark scripts
- Extended tutorial
- State-of-the-field comparison plan (JCat, OPTIMIZER, Codon Harmonizer, IDT) using documented, reproducible metrics only

---

## v3.8 — Stable Release

- All profiles have documented validation status; experimental profiles remain gated unless reviewed evidence supports activation
- Full documentation audit
- Final polish

---

## v4.0 — ML Engine *(data-conditional)*

No public ML engine is scheduled until sufficient, non-confidential validation data and benchmark evidence are available.

- ML engine exposed only after compatibility, cost, reproducibility, and claim-boundary review
- Training pipeline based on curated, non-confidential wet-lab validation summaries
- Rule-based vs ML benchmark comparison using documented public metrics
- No committed timeline — depends on wet-lab data accumulation

---

## Wet-lab Validation *(ongoing — not version-bound)*

Wet-lab validation is an ongoing, open-ended process independent of software versioning.
Public-safe summaries may be tracked as GitHub Issues under the [Wet-lab Validation](https://github.com/eijex/factorforge-cds/milestone/9) milestone. Private or sensitive feedback should not be submitted through public GitHub Issues.

- Submit only public-safe, non-confidential summaries via public [GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml). Send sensitive or private summaries by email instead.
- See [VALIDATION.md](VALIDATION.md) for submission format and current status

---

## Out of Scope

- Animal or clinical expression systems
- Protein structure prediction (use AlphaFold / ESMFold directly)
- CDS design for microbial hosts (E. coli, yeast)
