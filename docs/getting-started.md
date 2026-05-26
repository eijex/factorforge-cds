# Getting Started

## Requirements

- Python 3.10+

## Installation

```bash
pip install factorforge-cds
```

Experimental ML research modules (ESM2 + BART decoder):

```bash
pip install "factorforge-cds[ml]"
```

For development:

```bash
git clone https://github.com/eijex/factorforge-cds.git
cd factorforge
pip install -e ".[dev]"
```

## Docker (local web app)

Run the full web interface locally — no data leaves your machine:

```bash
docker pull ghcr.io/eijex/factorforge-cds:latest
docker run -p 8080:8080 ghcr.io/eijex/factorforge-cds:latest
```

Then open [http://localhost:8080](http://localhost:8080).

## Quick Start

**CLI:**

```bash
factorforge optimize my_protein.fasta -o output.fasta
```

**Python API:**

```python
from factorforge.engines.profile.pipeline import OptimizationPipeline

pipeline = OptimizationPipeline(profile="balanced")
result = pipeline.run("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG...")
print(result.sequence)   # optimized CDS
print(result.metadata)   # CAI, GC%, scan results, domestication edits
```

## Updating

**pip:**

```bash
pip install --upgrade factorforge-cds
```

**Docker:**

```bash
docker pull ghcr.io/eijex/factorforge-cds:latest
```

**Git clone:**

```bash
git pull origin main
pip install -e ".[dev]"
```

Check your installed version:

```bash
pip show factorforge-cds
# or
factorforge --version
```

Release notes: [CHANGELOG](changelog.md)
