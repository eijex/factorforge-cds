# Output Format

Each optimized sequence includes:

| Field | Description |
|-------|-------------|
| **Optimized CDS** | Synonymous codon replacements only — AA identity 100% |
| **CAI score** | Codon adaptation index for the selected host codon table |
| **GC content** | Global and first-region (first 90 nt) |
| **Scan report** | PolyA signals, CpG/TpA hotspots, homopolymers, rare codon runs, restriction sites |
| **Domestication report** | BsaI/BsmBI and custom restriction sites removed, edit count |
| **Construct ID** | Reproducible hash for tracking and reproducibility |

## FASTA Output

Profile-engine CLI `--output` FASTA uses self-describing metadata headers:

```
>target_optimized|engine=profile|profile=balanced|cai=0.810|gc=52.40|score=0.730
ATGAAACAGTTGGTCCTGGGG...
```

Profile output does not emit an `objective` field. `objective` is reserved for
the DP feasibility engine, where `objective=feasibility_best` identifies the
DP objective; profile-engine output is identified by `engine=profile` and
`profile=<profile>`.

## GenBank Output

GenBank format includes all metadata as FEATURES annotations, suitable for review in sequence editing tools (Benchling, SnapGene, etc.).

## Validator Contract

The design pipeline emits an output only after the hard-fail validator passes:

| Check | Requirement |
|-------|-------------|
| Amino acid identity | 100% |
| Internal stop codons | 0 |
| Invalid codons | 0 |
| Sequence length | Correct (input AA × 3) |

If any check fails, no output is produced.
