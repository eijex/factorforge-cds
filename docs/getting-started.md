# Getting Started

## Requirements

- Python 3.10+

## Installation

### Via Bioconda (Recommended for Bioinformatics)

```bash
conda create -n factorforge-review --override-channels -c conda-forge -c bioconda --strict-channel-priority python=3.11 factorforge-cds
conda activate factorforge-review
```

Bioconda officially supports Linux and macOS, not native Windows:
see the [Bioconda platform FAQ](https://bioconda.github.io/faqs.html).
Use WSL Ubuntu for this route on Windows. A local
installation does not require Neon, GCP, or a GPU for the DP optimization path.
Keep the installed package environment separate from a development checkout.

The Bioconda 3.5.11 package installed on native Windows in a limited test, but
does not provide a working Windows `factorforge.exe` launcher. The tested
workaround is `python -m factorforge` (for example,
`python -m factorforge optimize input.fasta -o output.fasta`). This observation
does not establish native Windows Bioconda support or validate optional research engines.

### Reviewer installed-package check

Bioconda 3.5.11 was smoke-tested on Ubuntu/WSL x86_64 with Python 3.11.
To replay that version, replace `factorforge-cds` in the installation command
with `factorforge-cds=3.5.11`. Linux/macOS are configured in the separate
`Bioconda installed package` CI workflow; configuration alone is not a passing
CI result or a claim that macOS has been tested locally.

After installing, run the maintained checker from a source download:

```bash
python -I /absolute/path/to/factorforge-cds/scripts/smoke_installed_distribution.py --output /tmp/factorforge-smoke
```

Do not run `pip install -e .` in that environment. The checker rejects a module
outside the environment or without a matching Bioconda installation record.
It runs the installed CLI in a temporary directory and records version/help,
default-versus-explicit DP output agreement, amino-acid identity and invalid-input
rejection. Its JSON report records the exact installed version and provenance,
including failures. It is a small installation smoke check, not the complete
software test suite or a biological benchmark; 3DP/4DP and live web behavior are
outside its scope. The CI checks the latest available Bioconda package, which
may lag behind the checkout or a newly published PyPI release.

### Via PyPI

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

## Sequence Length Limits

| Mode | Max length | Notes |
|------|------------|-------|
| Web API (protein input) | 5,000 aa | Vercel serverless timeout constraint |
| Web API (DNA/CDS input) | 15,000 bp | ≈ 5,000 aa |
| CLI / Docker | No limit | Local computation, install via `pip install factorforge-cds` |

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

## Eijex MCP

Use FactorForge from any MCP-compatible client:

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

Available tools:

| Tool | Description |
|------|-------------|
| `factorforge_cds_optimize` | Optimize a single protein sequence |
| `factorforge_cds_compare` | Compare multiple profiles side-by-side |
| `factorforge_cds_batch` | Optimize up to 20 sequences at once |

See [mcp.eijex.com](https://mcp.eijex.com) for the current public tool list.

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
