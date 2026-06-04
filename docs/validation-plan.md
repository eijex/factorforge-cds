# Validation Plan

FactorForge validation should cover implementation correctness, regression
stability, host-profile schema quality, score calibration, rule-audit behavior,
and documentation consistency.

## Validation Levels

| Level | What it checks | When |
|-------|---------------|------|
| Unit tests | Individual functions for scoring, rule engine, and reverse translator behavior | Every commit |
| Regression tests | Full optimization output for benchmark proteins | Every commit |
| Host profile schema | YAML schema validation for new profiles | On host profile changes |
| Scoring calibration snapshots | Composite score values for reference sequences | On scoring changes |
| Rule audit snapshots | Rule hit counts for benchmark sequences | On rule engine changes |
| Documentation consistency | Doc values match live code constants | On scoring or rule changes |

## Benchmark Protein Set

Regression and snapshot checks should use the benchmark protein set in:

```text
tests/fixtures/benchmark_proteins.py
```

The benchmark set is intended for stable sequence-design checks. It should be
safe, synthetic, and broad enough to exercise short sequences, reporter-like
sequences, antibody-like sequences, low-complexity regions, and cysteine-rich
regions.

## Snapshot Expectations

Scoring snapshots should record:

- Selected host and profile.
- CAI, GC, MFE status, and composite score.
- GC band boundaries and decay width.
- Active scoring weights after normalization.
- Optional scoring components that were unavailable or disabled.

Rule audit snapshots should record:

- Rule categories scanned.
- Thresholds and window sizes used.
- Finding counts by category.
- Representative context windows for deterministic checks.
- Whether repair was attempted and whether it succeeded.

## Documentation Consistency

Documentation consistency checks should compare published docs against live
constants for:

- Homopolymer expression and synthesis thresholds.
- GC band boundaries and decay width.
- Composite scoring weights.
- Rare codon run threshold and minimum run length.
- Implemented versus planned rule categories.
- Implemented versus planned MCP servers and agent skills.

Any scoring or rule change that can alter user-visible output should update
documentation in the same change set.
