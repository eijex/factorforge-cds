# FactorForge

**Open-source constraint-based CDS design engine for *Nicotiana benthamiana* expression workflows.**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-3.1.3-green.svg)](https://github.com/eijex/factorforge-cds/releases)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![Web App](https://img.shields.io/badge/web-factorforge--cds.vercel.app-brightgreen.svg)](https://factorforge-cds.vercel.app)

FactorForge optimizes protein sequences into *N. benthamiana*-compatible CDS by maximizing CAI, controlling GC content, eliminating PolyA signals, and producing MoClo/Golden Gate-ready constructs.

**→ [Full Documentation](https://eijex.github.io/factorforge-cds/)**

---

## Quick Start

```bash
pip install factorforge-cds
factorforge optimize my_protein.fasta -o output.fasta
```

Or use the **[web app](https://factorforge-cds.vercel.app)** — no installation required.

---

## Access Options

| Method | Description | Link |
|--------|-------------|------|
| **Web App** | No installation, demo & light use | [factorforge-cds.vercel.app](https://factorforge-cds.vercel.app) |
| **CLI / Python** | Local use, batch processing, data privacy | `pip install factorforge-cds` |
| **Docker** | Full web interface locally | `docker pull ghcr.io/eijex/factorforge-cds:latest` |

---

## ⚠️ Validation Status

FactorForge predictions are **in-silico only** and have not been experimentally validated in wet-lab conditions. See [Validation](https://eijex.github.io/factorforge-cds/validation/) and [VALIDATION.md](VALIDATION.md).

---

## Citing

```
FactorForge v3.1.3 (2026). Open-source constraint-based CDS design engine.
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
- **Web** — [factorforge-cds.vercel.app](https://factorforge-cds.vercel.app)
