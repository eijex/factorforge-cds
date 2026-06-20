# FactorForge

**Open-source constraint-based CDS design and pre-synthesis sequence review engine for plant CDS workflows, with primary support for *Nicotiana benthamiana* (Tobacco BY-2: experimental).**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](https://github.com/eijex/factorforge-cds/blob/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)
[![CI](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml/badge.svg)](https://github.com/eijex/factorforge-cds/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/eijex/factorforge-cds/branch/main/graph/badge.svg)](https://codecov.io/gh/eijex/factorforge-cds)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20407330.svg)](https://doi.org/10.5281/zenodo.20407330)

FactorForge performs profile-guided CDS design with CAI/GC metrics, PolyA-signal screening, and Golden Gate/MoClo-aware checks. It supports reproducible in-silico CDS candidate generation and pre-synthesis sequence review for plant CDS workflows. Primary support: *N. benthamiana* (agroinfiltration). Experimental host context: Tobacco BY-2 (`--host by2`).

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
| **[Eijex MCP](https://mcp.eijex.com)** | MCP-compatible agent access |

---

## Engineering Benchmark

Benchmarked on *N. benthamiana* SGN CDS (v3.2.3, balanced profile, N=49,257 sequences, seed=320):

| Metric | Value | Target |
|--------|-------|--------|
| CAI (mean) | 0.94 | ≥ 0.75 |
| GC in target range | 98.2% | 55–65% |
| AA identity | 100% | 100% |
| Validator pass rate | 100% | 100% |

> **Note**: The figures above are an *engineering reference benchmark* — a reproducible measurement against a specific version, profile, and dataset. A *formal comparative benchmark* against other CDS design tools is not yet available.

Reproduce the benchmark:

```bash
python scripts/benchmark.py --n 100  # quick validation (~30s)
python scripts/benchmark.py          # full benchmark (N. benthamiana SGN CDS)
```

Reproducible benchmark foundation: see `benchmarks/README.md`.

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
    FactorForge outputs are **in-silico only** and have not been experimentally validated in wet-lab conditions. They support reviewability and reproducibility, not guarantees of expression, yield, synthesis acceptance, folding, glycosylation, regulatory approval, or downstream biological performance. See [Validation](validation.md) for details.

## Paper-to-product journey

FactorForge's public roadmap keeps the paper/research-software track and product track connected but claim-bounded. The research track focuses on reproducible in-silico CDS candidate generation, documented benchmarks, and explicit limitations. The product track extends those artifacts toward a pre-synthesis review harness: assembly-relevant checks, host/profile maturity labels, quality-risk-aware annotations, and privacy-aware feedback records.

These extensions are future-facing review aids; they do not establish wet-lab performance, synthesis acceptance, biosecurity compliance, or regulatory approval.

---

## Share Wet-lab Results

Used FactorForge in the lab? Share public-safe feedback to help improve the tool:

- [Share Wet-lab Results (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) — public-safe coarse summaries only
- Email `eijex.lab@gmail.com` — private or sensitive summaries
