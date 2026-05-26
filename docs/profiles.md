# Optimization Profiles

Profiles apply to the profile engine (`--engine profile`). The `v2` engine key remains available as a compatibility alias.

| Profile | Description |
|---------|-------------|
| `balanced` | CAI + GC balance — default for most use cases |
| `high_cai` | Maximum codon adaptation index |
| `gc_target` | Target GC 42.5% for *N. benthamiana* |
| `viral_delivery` | Adjusted for TRV viral vector delivery |

## Usage

```bash
factorforge optimize input.fasta -e profile -p balanced -o output.fasta
factorforge optimize input.fasta -e profile -p high_cai -o output.fasta
factorforge optimize input.fasta -e profile -p viral_delivery -o output.fasta
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
