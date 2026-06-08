# FactorForge

**Open-source constraint-based CDS design engine for sequence-level CDS design, with primary support for *Nicotiana benthamiana* (Tobacco BY-2: experimental).**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407331.svg)](https://doi.org/10.5281/zenodo.20407331)
[![Web App](https://img.shields.io/badge/web-factorforge.eijex.com-brightgreen.svg)](https://factorforge.eijex.com)

FactorForge performs profile-guided CDS design with CAI/GC metrics, PolyA-signal screening, and Golden Gate/MoClo-aware checks. Primary support: *N. benthamiana* (agroinfiltration). Experimental host context: Tobacco BY-2 (`--host by2`).

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
| **Eijex MCP** | MCP-compatible agent access | [mcp.eijex.com](https://mcp.eijex.com) |

---

## Repository Structure

The supported production engine is the deterministic profile engine under:

```text
src/factorforge/engines/profile/
```

Historical implementation tracks are preserved under `archive/` for provenance
and are not imported by the installed package or exposed as supported engines.

---

## ⚠️ Validation Status

FactorForge outputs are **in-silico only** and have not been experimentally validated in wet-lab conditions. See [Validation](https://eijex.github.io/factorforge-cds/validation/) and [VALIDATION.md](VALIDATION.md).

---

## Citing

```
FactorForge v3.1.9 (2026). Open-source constraint-based CDS design engine.
Eijex. https://github.com/eijex/factorforge-cds
```

---

## Maintainer

Mun-Kyu Kim ([@eijex](https://github.com/eijex))

## License

GNU Affero General Public License v3.0 — see [LICENSE](LICENSE).

**Disclaimer:** FactorForge is provided for research purposes only. Outputs are computational and have not been experimentally validated.

---

## Get in Touch

- **Docs** — [eijex.github.io/factorforge-cds](https://eijex.github.io/factorforge-cds/)
- **Wet-lab Results** — Public-safe validation summaries are welcome. Do not submit raw sequences, confidential construct details, internal batch IDs, patient data, private contact information, exact process parameters, or confidential partner/customer data. See [VALIDATION.md](VALIDATION.md) before submitting.
- **GitHub Issues** — bugs, features: [github.com/eijex/factorforge-cds/issues](https://github.com/eijex/factorforge-cds/issues)
- **Email** — eijex.lab@gmail.com
- **FactorForge** — [factorforge.eijex.com](https://factorforge.eijex.com)
- **Lab** — [www.eijex.com](https://www.eijex.com)
