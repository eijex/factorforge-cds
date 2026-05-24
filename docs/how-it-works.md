# How It Works

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
   CpG/TpA dinucleotide hotspots,
   repeat sequences, rare codon runs,
   forbidden restriction sites
        │
        ▼
3. Domestication
   Removes Golden Gate / MoClo-incompatible
   BsaI / BsmBI recognition sites via silent edits
   Optional custom restriction sites removed
   CpG/TpA reduction with CAI-budgeted balanced mode
        │
        ▼
4. Output
   Optimized CDS — FASTA or GenBank
   with full metrics and scan report
```

## Engines

| Engine | Flag | Description |
|--------|------|-------------|
| `dp` | `--engine dp` | DP feasibility engine (default) — constraint-based, deterministic |
| `v2` | `--engine v2` | Rule-based optimizer — profile-driven CAI maximization |

## Design Objectives

| Objective | Description |
|-----------|-------------|
| `feasibility_best` | Best achievable CAI/GC under synonymous constraints (default) |
| `gc_target` | GC-focused rule-based output |
| `high_cai` | CAI-focused rule-based output |
