# FactorForge

**Open-source constraint-based CDS design engine for *Nicotiana benthamiana* expression workflows.**

[![License](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](https://github.com/eijex/factorforge-cds/blob/main/LICENSE)
[![Version](https://img.shields.io/badge/version-3.1.4-green.svg)](https://github.com/eijex/factorforge-cds/releases)
[![PyPI](https://img.shields.io/pypi/v/factorforge-cds.svg)](https://pypi.org/project/factorforge-cds/)

FactorForge optimizes protein sequences into *N. benthamiana*-compatible CDS by maximizing CAI, controlling GC content, eliminating PolyA signals, and producing MoClo/Golden Gate-ready constructs.

---

## Quick Start

```bash
pip install factorforge-cds
factorforge optimize my_protein.fasta -o output.fasta
```

Or use the **[web app](https://factorforge-cds.vercel.app)** — no installation required.

---

## Access Options

| Method | Description |
|--------|-------------|
| **[Web App](https://factorforge-cds.vercel.app)** | No installation, demo & light use |
| **CLI / Python** | Local use, batch processing, data privacy |
| **Docker** | Full web interface locally, no data leaves your machine |

---

## Performance

Benchmarked on *N. benthamiana* SGN CDS (v3.1.4, balanced profile, N=3,876 sequences):

| Metric | Value | Target |
|--------|-------|--------|
| CAI (mean) | 0.76 | ≥ 0.75 |
| GC% (mean) | 59.77% | 55–65% |
| AA identity | 100% | 100% |
| Validator pass rate | 100% | 100% |

---

## Supported Hosts

| Host | Status |
|------|--------|
| *Nicotiana benthamiana* | ✅ Supported |
| *Wolffia globosa* | 🔶 Codon table available, coming soon |
| Other plant hosts | 📋 Planned |

---

!!! warning "Validation Status"
    FactorForge predictions are **in-silico only** and have not been experimentally validated in wet-lab conditions. See [Validation](validation.md) for details.
