# FactorForge

**Open-source constraint-based CDS design and pre-synthesis sequence review for plant CDS workflows, with primary support for *Nicotiana benthamiana* (Tobacco BY-2: experimental).**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407330.svg)](https://doi.org/10.5281/zenodo.20407330)
[![Web App](https://img.shields.io/badge/web-factorforge.eijex.com-brightgreen.svg)](https://factorforge.eijex.com)
[![Buy Me A Coffee](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-Donate-yellow?style=flat&logo=buy-me-a-coffee)](https://www.buymeacoffee.com/eijex)


FactorForge performs deterministic CDS design with CAI/GC metrics, PolyA-signal screening, and Golden Gate/MoClo-aware checks. It is positioned as a pre-synthesis review harness: it helps teams generate reproducible CDS candidates, inspect assembly-relevant sequence constraints, and package design metadata before downstream synthesis, cloning, or experimental review. Primary support: *N. benthamiana* (agroinfiltration). Experimental host context: Tobacco BY-2 (`--host by2`).

The v3.5.0 release-candidate line uses independently versioned engines:

| Generation | Engine | Version | Availability |
|---|---|---:|---|
| Gen 1 | Rule/profile | 1.0.0 | Stable, public |
| Gen 2 | DP v2 | 2.0.0 | Release candidate, public |
| Gen 3 | sLLM Hybrid | 0.1.0-preview.1 | Feature-gated research preview; no trained model is claimed |

**→ [Full Documentation](https://eijex.github.io/factorforge-cds/)** · **[Roadmap](https://eijex.github.io/factorforge-cds/roadmap/)**

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

The supported deterministic engines are the profile engine and DP v2 under:

```text
src/factorforge/engines/profile/
src/factorforge/engines/dp_v2.py
```

Historical implementation tracks are preserved under `archive/` for provenance
and are not imported by the installed package or exposed as supported engines.

---

## ⚠️ Validation Status

FactorForge outputs are **in-silico only** and have not been experimentally validated in wet-lab conditions. These checks support reviewability and reproducibility; they do not guarantee expression, yield, synthesis acceptance, folding, glycosylation, regulatory approval, or downstream biological performance. See [Validation](https://eijex.github.io/factorforge-cds/validation/) and [VALIDATION.md](VALIDATION.md).

---

## Citing

```
FactorForge v3.4.5 (2026). Open-source constraint-based CDS design and sequence review.
Eijex. https://github.com/eijex/factorforge-cds
```

---

## Maintainer

Mun-Kyu Kim ([@eijex](https://github.com/eijex))

## Optional PostgreSQL support

Install `pip install "factorforge-cds[postgres]"` only when using the PostgreSQL
connector. Explicit local SQLite checkpoints do not require a PostgreSQL driver.
Database integration tests require an explicitly configured `DATABASE_URL`.

## License

GNU Affero General Public License v3.0 — see [LICENSE](LICENSE).

**Disclaimer:** FactorForge is provided for research purposes only. Outputs are computational and have not been experimentally validated.

---

## Get in Touch

- **Sponsor** — Support our research via [GitHub Sponsors](https://github.com/sponsors/eijex) or [PayPal](https://paypal.me/munkyukim)
- **Docs** — [eijex.github.io/factorforge-cds](https://eijex.github.io/factorforge-cds/)
- **Wet-lab Feedback** — Public-safe feedback summaries are welcome via [Share Wet-lab Feedback (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml). Do not submit raw sequences, confidential construct details, internal batch IDs, patient data, private contact information, exact process parameters, or confidential partner/customer data. Email `eijex.lab@gmail.com` for private or sensitive summaries. See [VALIDATION.md](VALIDATION.md) before submitting.
- **GitHub Issues** — bugs, features: [github.com/eijex/factorforge-cds/issues](https://github.com/eijex/factorforge-cds/issues)
- **Email** — eijex.lab@gmail.com
- **FactorForge** — [factorforge.eijex.com](https://factorforge.eijex.com)
- **Eijex MCP** — [mcp.eijex.com](https://mcp.eijex.com)
- **Lab** — [www.eijex.com](https://www.eijex.com)

