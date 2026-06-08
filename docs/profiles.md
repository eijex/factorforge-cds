# Optimization Profiles

Profiles apply to the profile engine (`--engine profile`).

## Supported Hosts

| Host flag | Species | Codon table | Notes |
|-----------|---------|-------------|-------|
| `nbenthamiana` (default) | *Nicotiana benthamiana* | 3,876 SGN CDS | Stable sequence-level design profile for agroinfiltration workflows |
| `by2` ⚠️ Experimental | *N. tabacum* BY-2 | 1,534 Kazusa CDS (2007) | **Experimental.** Uses *N. tabacum* codon usage as proxy. Not wet-lab validated for BY-2 expression. |

!!! warning "BY-2 host is experimental"
    The `by2` host profile uses *N. tabacum* codon usage data as a proxy for Tobacco BY-2 suspension cells. It has not been wet-lab validated for BY-2 expression performance. Use for exploratory design only.

## Stable Profiles

These four profiles are fully supported and available via CLI, Python API, web app, and MCP.

| Profile | Description |
|---------|-------------|
| `balanced` | CAI + GC balance — default for most sequence-design reviews |
| `high_cai` | CAI-focused synonymous CDS candidate |
| `gc_target` | Targets a configurable GC percentage, defaulting to the host-profile midpoint (60% for *N. benthamiana*). Pass an explicit target to drive GC higher or lower. |
| `assembly_friendly` | Golden Gate / MoClo workflows — avoids BsaI/BpiI Type IIS restriction sites via synonymous substitution. Does not yet score local GC uniformity or repeat patterns. |

## Usage

```bash
factorforge optimize input.fasta -e profile -p balanced -o output.fasta
factorforge optimize input.fasta -e profile -p high_cai -o output.fasta
factorforge optimize input.fasta -e profile -p gc_target -o output.fasta
factorforge optimize input.fasta -e profile -p assembly_friendly -o output.fasta
factorforge optimize input.fasta -e profile -p balanced --host by2 -o output.fasta
```

## Python API

```python
from factorforge.engines.profile.pipeline import OptimizationPipeline

pipeline = OptimizationPipeline(profile="high_cai")
result = pipeline.run("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG...")
```

## Profile Selection Guide

- **Default *N. benthamiana* CDS design** → `balanced`
- **CAI-focused comparison** → `high_cai`
- **GC-sensitive downstream processes** → `gc_target` (pass an explicit GC target if you need a value other than the host midpoint)
- **Golden Gate / MoClo assembly** → `assembly_friendly`

## Not Public API

The following profiles are **not part of the public API or MCP interface** and are not supported for general use. They may change or be removed without notice.

| Profile | Status | Notes |
|---------|--------|-------|
| `viral_delivery` | Not public | Adjusted for TRV viral vector delivery. Not validated. |
| `ramp` | Not public | 5′ ramp design heuristic. Under development. |
