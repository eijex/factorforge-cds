# Benchmark Suite Report: Engineering Benchmark Suite
**Run Timestamp**: 2026-08-29T14:35:06.721958

## Summary (Pass Rate)
*(Format: Passed / Total Targets)*
- **profile**: 0/2
- **dp**: 0/2
- **lm**: 0/2

## Average Metrics
| Engine | Average CAI | Average GC% |
|--------|-------------|-------------|
| profile | 0.931 | 33.9% |
| dp | 0.940 | 39.9% |
| lm | 0.875 | 39.5% |

## Target Results
### Target: eGFP | Engine: profile | Passed: ❌
- **Runtime**: 7.17s
- **GC%**: 31.7%
- **CAI**: 0.944

**Hard/Gate Failures**
- gc_content: GC 31.67% is outside target 40.0%-47.0%

**Warnings**
- polya_signal: Found 4 PolyA signals
- homopolymer: Found 3 homopolymer runs

**Passed Invariants**
- frame_valid
- aa_identity
- internal_stop
- type_iis_bsai
- type_iis_bpii
- type_iis_bsmbi
- direct_repeat

### Target: eGFP | Engine: dp | Passed: ❌
- **Runtime**: 0.14s
- **GC%**: 39.9%
- **CAI**: 0.916

**Hard/Gate Failures**
- gc_content: GC 39.86% is outside target 40.0%-47.0%

**Passed Invariants**
- frame_valid
- aa_identity
- internal_stop
- type_iis_bsai
- type_iis_bpii
- type_iis_bsmbi
- polya_signal
- homopolymer
- direct_repeat

### Target: eGFP | Engine: lm | Passed: ❌
- **Runtime**: 4.08s
- **GC%**: 39.3%
- **CAI**: 0.854

**Hard/Gate Failures**
- gc_content: GC 39.31% is outside target 40.0%-47.0%

**Warnings**
- polya_signal: Found 1 PolyA signals

**Passed Invariants**
- frame_valid
- aa_identity
- internal_stop
- type_iis_bsai
- type_iis_bpii
- type_iis_bsmbi
- homopolymer
- direct_repeat

### Target: Humira_HC | Engine: profile | Passed: ❌
- **Runtime**: 0.03s
- **GC%**: 36.2%
- **CAI**: 0.917

**Hard/Gate Failures**
- gc_content: GC 36.21% is outside target 40.0%-47.0%

**Warnings**
- polya_signal: Found 3 PolyA signals
- homopolymer: Found 5 homopolymer runs

**Passed Invariants**
- frame_valid
- aa_identity
- internal_stop
- type_iis_bsai
- type_iis_bpii
- type_iis_bsmbi
- direct_repeat

### Target: Humira_HC | Engine: dp | Passed: ❌
- **Runtime**: 1.15s
- **GC%**: 40.0%
- **CAI**: 0.964

**Hard/Gate Failures**
- gc_content: GC 39.97% is outside target 40.0%-47.0%
- type_iis_bsai: Found forbidden BsaI sites: ['GGTCTC']
- type_iis_bpii: Found forbidden BpiI sites: ['GAAGAC', 'GTCTTC']

**Warnings**
- homopolymer: Found 2 homopolymer runs

**Passed Invariants**
- frame_valid
- aa_identity
- internal_stop
- type_iis_bsmbi
- polya_signal
- direct_repeat

### Target: Humira_HC | Engine: lm | Passed: ❌
- **Runtime**: 0.17s
- **GC%**: 39.7%
- **CAI**: 0.896

**Hard/Gate Failures**
- gc_content: GC 39.75% is outside target 40.0%-47.0%

**Warnings**
- polya_signal: Found 3 PolyA signals

**Passed Invariants**
- frame_valid
- aa_identity
- internal_stop
- type_iis_bsai
- type_iis_bpii
- type_iis_bsmbi
- homopolymer
- direct_repeat
