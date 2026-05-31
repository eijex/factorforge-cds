# Getting Started

## Requirements

- Python 3.10+

## Installation

```bash
pip install factorforge-cds
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

## Eijex MCP (AI agent access)

Use FactorForge directly from Claude Code, Cursor, or any MCP-compatible client:

```json
{
  "mcpServers": {
    "eijex": {
      "type": "http",
      "url": "https://mcp.eijex.com/api/mcp"
    }
  }
}
```

Available tools: `factorforge_cds_optimize`, `factorforge_cds_compare`, `factorforge_cds_batch`

See [mcp.eijex.com](https://mcp.eijex.com) for the full tool list.

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
