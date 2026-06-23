# How It Works

FactorForge runs a deterministic constraint-based pipeline in four stages:

```
Protein sequence (FASTA or plain text)
        │
        ▼
1. Reverse Translation
   Selects synonymous codons to favor CAI
   against the selected host codon usage table
   (N. benthamiana by default; --host by2 for Tobacco BY-2)
        │
        ▼
2. Rule Scan (9 advisory scanners, run by default)
   Detects PolyA signals, AU-rich elements,
   AT-rich runs, homopolymers, tandem repeats,
   local GC extremes, splice-like motifs,
   CpG/TpA dinucleotide hotspots, rare codon runs
   (advisory findings only — never gating)
        │
        ▼
3. Domestication
   Removes Golden Gate / MoClo-incompatible
   BsaI / BsmBI recognition sites via silent edits
   (halts if unresolvable)
   Optional custom restriction sites removed
   CpG/TpA reduction with CAI-budgeted balanced mode
   MoClo overhang validity/collision check available
   when building a construct (opt-in, non-gating)
        │
        ▼
4. Output
   Designed CDS — FASTA or GenBank
   with full metrics and scan report
```

## Engines

| Engine | Flag | Description |
|--------|------|-------------|
| `dp` | `--engine dp` | DP feasibility engine (default) — constraint-based, deterministic |
| `profile` | `--engine profile` | Rule/profile optimizer — profile-driven CAI and constraint handling |

## Design Objectives

The DP engine (`--engine dp`) supports one `--objective`:

| Objective | Description |
|-----------|-------------|
| `feasibility_best` | Best achievable CAI/GC under synonymous constraints (default) |

The profile engine (`--engine profile`) supports these `--profile` values
(not DP `--objective` values — passing `gc_target` or `high_cai` to
`--objective` always fails, since the DP engine only implements
`feasibility_best`; see [Optimization Profiles](profiles.md)):

| Profile | Description |
|---------|-------------|
| `gc_target` | GC-focused rule-based output |
| `high_cai` | CAI-focused rule-based output |
