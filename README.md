# FactorForge

**Open-source constraint-based CDS design engine for *Nicotiana benthamiana* expression workflows.**

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-3.0.0-green.svg)](https://github.com/eijex/factorforge-cds/releases)
[![Web App](https://img.shields.io/badge/web-factorforge-cds.vercel.app-brightgreen.svg)](https://factorforge-cds.vercel.app)

FactorForge optimizes protein sequences into *N. benthamiana*-compatible CDS by maximizing CAI, controlling GC content, eliminating PolyA signals, and producing MoClo/Golden Gate-ready constructs.

---

## Quick Start

```bash
pip install factorforge-cds
factorforge optimize my_protein.fasta -o output.fasta
```

Or with Python:

```python
from factorforge.engines.v2.pipeline import OptimizationPipeline

pipeline = OptimizationPipeline(profile="balanced")
result = pipeline.run("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG...")
print(result.sequence)   # optimized CDS
print(result.metadata)   # CAI, GC%, scan results, domestication edits
```

---

## Access Options

| Method | Description | Link |
|--------|-------------|------|
| **Web App** | No installation, demo & light use | [factorforge-cds.vercel.app](https://factorforge-cds.vercel.app) |
| **CLI / Python** | Local use, batch processing, data privacy | `pip install factorforge-cds` |
| **Notebooks** | Training & experimentation on Colab / Kaggle | See [notebooks/](notebooks/) |

---

## How It Works

FactorForge runs a deterministic constraint-based pipeline in four stages:

```
Protein sequence (FASTA or plain text)
        │
        ▼
1. Reverse Translation
   Selects synonymous codons to maximize CAI
   against the N. benthamiana codon usage table
        │
        ▼
2. Rule Scan
   Detects PolyA signals, homopolymers,
   repeat sequences, rare codon runs,
   forbidden restriction sites
        │
        ▼
3. Domestication
   Removes Golden Gate / MoClo-incompatible
   BsaI / BsmBI recognition sites via silent edits
   Optional custom restriction sites can be removed
   by synonymous substitution when feasible
        │
        ▼
4. Output
   Optimized CDS — FASTA or GenBank
   with full metrics and scan report
```

### Optimization Profiles

| Profile | Description |
|---------|-------------|
| `balanced` | CAI + GC balance (default) |
| `high_cai` | Maximum codon adaptation |
| `gc_target` | Target GC 42.5% for N. benthamiana |
| `viral_delivery` | Adjusted for viral vector delivery |

---

## Performance

Benchmarked on *N. benthamiana* codon usage table (v2 engine, 3,876 sequences):

| Metric | Value | Target |
|--------|-------|--------|
| CAI (mean) | 0.80 | ≥ 0.75 |
| GC% (mean) | 42.54% | 40–55% |
| GC% (range) | 40.36–53.81% | 40–55% |
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

## Installation

**Requirements:** Python 3.10+

```bash
pip install factorforge-cds
```

Experimental ML research modules are available separately:

```bash
pip install "factorforge-cds[ml]"
```

> These modules (ESM2 + BART decoder) are not part of the stable v3.0.0 default optimizer.
> The default v3.0.0 engine is the constraint-based DP feasibility engine.

For development:

```bash
git clone https://github.com/eijex/factorforge-cds.git
cd factorforge
pip install -e ".[dev]"
```

### Docker (local web app)

Run the full web interface locally — no data leaves your machine:

```bash
docker pull ghcr.io/eijex/factorforge-cds:latest
docker run -p 8080:8080 ghcr.io/eijex/factorforge-cds:latest
```

Then open [http://localhost:8080](http://localhost:8080).

Or build from source:

```bash
git clone https://github.com/eijex/factorforge-cds.git
cd factorforge
docker build -t factorforge-cds .
docker run -p 8080:8080 factorforge-cds
```

---

## CLI Reference

```bash
# Basic optimization (DP feasibility engine, default)
factorforge optimize input.fasta -o output.fasta

# Rule-based engine with profile
factorforge optimize input.fasta -e v2 -p balanced -o output.fasta

# With MoClo construct template, GenBank output
factorforge optimize input.fasta -e v2 -p balanced \
  --template standard_expression -o output.gb --format genbank

# Custom GC target range
factorforge optimize input.fasta --gc-min 40 --gc-max 50 -o output.fasta

# List available engines
factorforge list-engines
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--engine`, `-e` | `dp` | Engine: `dp` (feasibility) or `v2` (rule-based) |
| `--profile`, `-p` | `balanced` | Optimization profile |
| `--objective` | `feasibility_best` | DP objective |
| `--gc-min` / `--gc-max` | 40 / 55 | GC% target range |
| `--format` | `fasta` or `genbank` | Output format |
| `--scan-mode` | `full` | Rule scan: `full` or `fast` |

---

## Output

Each optimized sequence includes:

- **Optimized CDS** — synonymous codon replacements only, AA identity 100%
- **CAI score** — codon adaptation index for *N. benthamiana*
- **GC content** — global and first-region
- **Scan report** — PolyA signals detected/fixed, homopolymers, rare codon runs, restriction sites
- **Domestication report** — BsaI/BsmBI and optional custom restriction sites removed, edit count
- **Construct ID** — reproducible hash for tracking

---

## ⚠️ Validation Status

FactorForge predictions are **in-silico only** and have not been experimentally
validated in wet-lab conditions.

We are actively seeking researchers to test these predictions. If you use
FactorForge in your experiments, we'd love to hear from you:

- Did the optimized sequence express well?
- How did CAI / GC% correlate with actual expression levels?
- Any unexpected results?

**Share your results** → [GitHub Issues](https://github.com/eijex/factorforge-cds/issues)
or email: eijex.lab@gmail.com

Validated results will be credited in [VALIDATION.md](VALIDATION.md) and future releases.

---

## 🛠️ Developed With

This project was built using the following tools and platforms:

| Tool | Role |
|------|------|
| [Claude](https://claude.ai) / [Claude Code](https://claude.ai/code) (Anthropic) | Architecture design, domain analysis, code review |
| [Codex](https://github.com/openai/codex) (OpenAI) | Code generation and implementation |
| [Google Colab](https://colab.research.google.com) | ML training experiments |
| [Kaggle](https://www.kaggle.com) | ML training experiments |
| [ESM2](https://github.com/facebookresearch/esm) (Meta) | Protein language model (encoder) |
| [PyTorch](https://pytorch.org) | ML framework |
| [Conda](https://docs.conda.io) / [Miniconda](https://docs.anaconda.com/miniconda/) | Environment management |
| [Vercel](https://vercel.com) | Web deployment |
| [GitHub](https://github.com) | Version control and open-source hosting |
| [HuggingFace](https://huggingface.co) | Model ecosystem |
| [BioPython](https://biopython.org) | Biological sequence processing |

---

## Citing

If you use FactorForge in your research, please cite:

```
FactorForge v3.0.0 (2026). Open-source constraint-based CDS design engine.
Eijex. https://github.com/eijex/factorforge-cds
```

*A citable publication is in preparation. Until then, please cite the GitHub repository.*

---

## License & Disclaimer

FactorForge source code is licensed under the [Apache License 2.0](LICENSE).

**Disclaimer:** FactorForge is provided for research purposes only. Predictions
are computational and have not been experimentally validated. The authors make
no warranties regarding expression outcomes in wet-lab settings. Use at your
own discretion.

---

## Get in Touch

- **GitHub Issues** — bug reports, feature requests, wet-lab results: [github.com/eijex/factorforge-cds/issues](https://github.com/eijex/factorforge-cds/issues)
- **Email** — collaborations, feedback, questions: eijex.lab@gmail.com
- **Web** — [factorforge-cds.vercel.app](https://factorforge-cds.vercel.app)
