# Optimization Profiles

Profiles apply to the v2 rule-based engine (`--engine v2`).

| Profile | Description |
|---------|-------------|
| `balanced` | CAI + GC balance — default for most use cases |
| `high_cai` | Maximum codon adaptation index |
| `gc_target` | Target GC 42.5% for *N. benthamiana* |
| `viral_delivery` | Adjusted for TRV viral vector delivery |

## Usage

```bash
factorforge optimize input.fasta -e v2 -p balanced -o output.fasta
factorforge optimize input.fasta -e v2 -p high_cai -o output.fasta
factorforge optimize input.fasta -e v2 -p viral_delivery -o output.fasta
```

## Python API

```python
from factorforge.engines.v2.pipeline import OptimizationPipeline

pipeline = OptimizationPipeline(profile="high_cai")
result = pipeline.run("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG...")
```

## Profile Selection Guide

- **General expression** → `balanced`
- **Maximizing expression level** → `high_cai`
- **GC-sensitive downstream processes** → `gc_target`
- **Viral vector (TRV, TMV)** → `viral_delivery`
