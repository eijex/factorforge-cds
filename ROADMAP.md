# Roadmap

FactorForge development direction and planned work.

> Current release: **v3.1.4** — [Changelog](CHANGELOG.md) · [Releases](https://github.com/eijex/factorforge-cds/releases)

---

## Current Focus (v3.1.x)

Wet-lab validation and stability. The core engine is feature-complete; this track collects experimental results and applies targeted fixes based on real expression data.

- Wet-lab feedback integration (Google Form → GitHub Issues)
- Bug fixes and metric corrections based on community reports
- Documentation improvements

---

## Planned (v3.2)

| Feature | Description |
|---------|-------------|
| **Wolffia globosa support** | Codon table available; enable optimization for this host |
| **Additional plant hosts** | Codon tables and profile tuning for other expression systems |
| **Batch optimization API** | Multi-sequence input for pipeline use cases |
| **Extended DP engine** | Expose more DP feasibility parameters via CLI and API |

---

## ML Research Track (v4)

A future ML-based CDS design engine is planned once sufficient wet-lab validation data is available to train and evaluate against. The v3-alpha ML prototype (archived under `archive/v3-ml-prototype/`) informs the architecture for this track.

This track has no committed timeline — it depends on wet-lab data accumulation.

---

## Wet-lab Validation

Wet-lab validation is an ongoing, open-ended process independent of software versioning.

- Submit results via [Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform?usp=header) or [GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml)
- See [VALIDATION.md](VALIDATION.md) for submission format and current status

---

## Out of Scope

- Animal or clinical expression systems
- Codon optimization for non-plant hosts
- Protein structure prediction (use AlphaFold / ESMFold directly)
