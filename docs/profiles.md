# Optimization Profiles

Profiles apply to the profile engine (`--engine profile`).

## Supported Hosts

| Host flag | Species | Codon table | Notes |
|-----------|---------|-------------|-------|
| `nbenthamiana` (default) | *Nicotiana benthamiana* | 3,876 SGN CDS | Stable. Recommended for agroinfiltration |
| `by2` ⚠️ Experimental | *N. tabacum* BY-2 | 1,534 Kazusa CDS (2007) | **Experimental.** Uses *N. tabacum* codon usage as proxy. Not wet-lab validated for BY-2 expression. |

!!! warning "BY-2 host is experimental"
    The `by2` host profile uses *N. tabacum* codon usage data as a proxy for Tobacco BY-2 suspension cells. It has not been wet-lab validated for BY-2 expression performance. Use for exploratory design only.

| Profile | Description |
|---------|-------------|
| `balanced` | CAI + GC balance — default for most use cases |
| `high_cai` | Maximum codon adaptation index |
| `gc_target` | Target GC 42.5% for *N. benthamiana* |
| `viral_delivery` | Adjusted for TRV viral vector delivery |
| `ml_enhanced` ⚠️ Experimental | SynCodonLM-augmented scoring (`pip install factorforge-cds[ml]` required) |

## Usage

```bash
factorforge optimize input.fasta -e profile -p balanced -o output.fasta
factorforge optimize input.fasta -e profile -p high_cai -o output.fasta
factorforge optimize input.fasta -e profile -p viral_delivery -o output.fasta
factorforge optimize input.fasta -e profile -p balanced --host by2 -o output.fasta
```

## Python API

```python
from factorforge.engines.profile.pipeline import OptimizationPipeline

pipeline = OptimizationPipeline(profile="high_cai")
result = pipeline.run("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG...")
```

## Profile Selection Guide

- **General expression** → `balanced`
- **Maximizing CAI** → `high_cai`
- **GC-sensitive downstream processes** → `gc_target`
- **Viral vector (TRV, TMV)** → `viral_delivery`
- **ML-augmented scoring (experimental)** → `ml_enhanced` (requires `pip install factorforge-cds[ml]`)
