# Output Format

Each optimized sequence includes:

| Field | Description |
|-------|-------------|
| **Optimized CDS** | Synonymous codon replacements only — AA identity 100% |
| **CAI score** | Codon adaptation index for the selected expression host |
| **GC content** | Global and first-region (first 90 nt) |
| **Scan report** | PolyA signals, CpG/TpA hotspots, homopolymers, rare codon runs, restriction sites |
| **Domestication report** | BsaI/BsmBI and custom restriction sites removed, edit count |
| **Construct ID** | Reproducible hash for tracking and reproducibility |

## FASTA Output

```
>CD47_optimized | CAI=0.81 | GC=52.4% | construct_id=a3f7b2
ATGAAACAGTTGGTCCTGGGG...
```

## GenBank Output

GenBank format includes all metadata as FEATURES annotations, suitable for direct import into sequence editing tools (Benchling, SnapGene, etc.).

## Validator Contract

All outputs are guaranteed to pass the hard-fail validator:

| Check | Requirement |
|-------|-------------|
| Amino acid identity | 100% |
| Internal stop codons | 0 |
| Invalid codons | 0 |
| Sequence length | Correct (input AA × 3) |

If any check fails, no output is produced.
