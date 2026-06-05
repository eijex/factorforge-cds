# FactorForge

**Open-source constraint-based CDS design engine for plant expression workflows, with initial focus on *Nicotiana benthamiana* and Tobacco BY-2.**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](https://github.com/eijex/factorforge-cds/blob/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407331.svg)](https://doi.org/10.5281/zenodo.20407331)

FactorForge optimizes protein sequences into host-compatible CDS by maximizing CAI, controlling GC content, eliminating PolyA signals, and producing MoClo/Golden Gate-ready constructs. Supports *N. benthamiana* (agroinfiltration) and Tobacco BY-2 (`--host by2`, bioreactor/cGMP).

---

## Quick Start

```bash
pip install factorforge-cds
factorforge optimize my_protein.fasta -o output.fasta
```

Or use the **[web app](https://factorforge.eijex.com)** — no installation required.

---

## Access Options

| Method | Description |
|--------|-------------|
| **[Web App](https://factorforge.eijex.com)** | No installation, demo & light use |
| **CLI / Python** | Local use, batch processing, data privacy |
| **Docker** | Full web interface locally, no data leaves your machine (`docker pull ghcr.io/eijex/factorforge-cds:latest`) |
| **[Eijex MCP](https://mcp.eijex.com)** | AI agent access via Claude Code, Cursor, and other MCP clients |

---

## Performance

Benchmarked on *N. benthamiana* SGN CDS (v3.1.4, balanced profile, N=3,876 sequences):

| Metric | Value | Target |
|--------|-------|--------|
| CAI (mean) | 0.76 | ≥ 0.75 |
| GC% (mean) | 59.77% | 55–65% |
| AA identity | 100% | 100% |
| Validator pass rate | 100% | 100% |

Reproduce the benchmark:

```bash
python scripts/benchmark.py --n 100  # quick validation (~30s)
python scripts/benchmark.py          # full benchmark (N. benthamiana SGN CDS)
```

---

## Supported Hosts

| Host | Status |
|------|--------|
| *Nicotiana benthamiana* | ✅ Supported (`--host nbenthamiana`, default) |
| Tobacco BY-2 (*N. tabacum*) | ⚠️ Experimental (`--host by2`) |
| *Wolffia globosa* | 🔶 Codon table available, coming soon |
| Other plant hosts | 📋 Planned |

---

!!! warning "Validation Status"
    FactorForge predictions are **in-silico only** and have not been experimentally validated in wet-lab conditions. See [Validation](validation.md) for details.
