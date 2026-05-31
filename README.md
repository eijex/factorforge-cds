# FactorForge

**Open-source constraint-based CDS design engine for plant expression workflows, with initial focus on *Nicotiana benthamiana* and Tobacco BY-2.**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407331.svg)](https://doi.org/10.5281/zenodo.20407331)
[![Web App](https://img.shields.io/badge/web-factorforge.eijex.com-brightgreen.svg)](https://factorforge.eijex.com)

FactorForge optimizes protein sequences into host-compatible CDS by maximizing CAI, controlling GC content, eliminating PolyA signals, and producing MoClo/Golden Gate-ready constructs. Supports *N. benthamiana* (agroinfiltration) and Tobacco BY-2 (`--host by2`, bioreactor/cGMP workflows).

**→ [Full Documentation](https://eijex.github.io/factorforge-cds/)**

---

## Quick Start

```bash
pip install factorforge-cds
factorforge optimize my_protein.fasta -o output.fasta
```

Or use the **[web app](https://factorforge.eijex.com)** — no installation required.

---

## Access Options

| Method | Description | Link |
|--------|-------------|------|
| **Web App** | No installation, demo & light use | [factorforge.eijex.com](https://factorforge.eijex.com) |
| **CLI / Python** | Local use, batch processing, data privacy | `pip install factorforge-cds` |
| **Docker** | Full web interface locally | `docker pull ghcr.io/eijex/factorforge-cds:latest` |

---

## Repository Structure

The supported production engine is the deterministic profile engine under:

```text
src/factorforge/engines/profile/
```

Historical implementation tracks are preserved under `archive/` for provenance
and are not imported by the installed package or exposed as supported engines.

---

## Development History

FactorForge has gone through several implementation generations before the current public release:

| Generation | Status | Description |
|-----------|--------|-------------|
| **v1** — NBent_OptiCodon | Internal | Thesis-derived codon optimization baseline for *N. benthamiana* |
| **v2** — Rule-Based Engine | Internal → Production | Deterministic, constraint-aware design engine; became the foundation for the public release |
| **v3-alpha** — ML Prototype | Archived | ML-based design attempt; performance was insufficient for production use; preserved under `archive/v3-ml-prototype/` |
| **v3.0+** — Current release | Public | Open-source release of the matured v2 engine under `factorforge.engines.profile` |
| **v3.7+** — ML Engine | Planned | ML-based design as `--engine ml`; added once sufficient wet-lab data is available |

The `archive/` directory preserves all three earlier tracks for provenance. None are installed or exposed by the current package.

---

## ⚠️ Validation Status

FactorForge predictions are **in-silico only** and have not been experimentally validated in wet-lab conditions. See [Validation](https://eijex.github.io/factorforge-cds/validation/) and [VALIDATION.md](VALIDATION.md).

---

## Citing

```
FactorForge v3.1.6 (2026). Open-source constraint-based CDS design engine.
Eijex. https://github.com/eijex/factorforge-cds
```

*A citable publication is in preparation.*

---

## Contributors

| | Name | Role |
|--|------|------|
| 👤 | Mun-Kyu Kim ([@eijex](https://github.com/eijex)) | Author & maintainer |
| 🤖 | Claude (Anthropic) | Design, analysis, planning |
| 🤖 | Codex (OpenAI) | Implementation |

## License

GNU Affero General Public License v3.0 — see [LICENSE](LICENSE).

**Disclaimer:** FactorForge is provided for research purposes only. Predictions are computational and have not been experimentally validated.

---

## Get in Touch

- **Docs** — [eijex.github.io/factorforge-cds](https://eijex.github.io/factorforge-cds/)
- **Wet-lab Results** — [Submit via Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform?usp=header) (recommended) or [GitHub Issue](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml)
- **GitHub Issues** — bugs, features: [github.com/eijex/factorforge-cds/issues](https://github.com/eijex/factorforge-cds/issues)
- **Email** — eijex.lab@gmail.com
- **Web** — [factorforge.eijex.com](https://factorforge.eijex.com)
