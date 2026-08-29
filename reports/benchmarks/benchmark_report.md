# Benchmark Suite Report: Engineering Benchmark Suite
**Run Timestamp**: 2026-08-29T10:33:47.694658

## Summary (Pass Rate)
*(Format: Passed / Total Targets)*
- **profile**: 0/2
- **dp**: 0/2
- **lm**: 0/2

## Average Metrics
| Engine | Average CAI | Average GC% |
|--------|-------------|-------------|
| profile | 0.000 | 3404.4% |
| dp | 0.000 | 3322.7% |
| lm | 0.000 | 3952.7% |

## Target Results
### Target: eGFP | Engine: profile | Passed: ❌
- **Runtime**: 2.47s
- **GC%**: 3179.9%
- **Failing Checks**:
  - frame_valid (CheckEnforcement.HARD_FAIL): Sequence length (717) is not divisible by 3
  - aa_identity (CheckEnforcement.HARD_FAIL): AA Identity is 1.000
  - internal_stop (CheckEnforcement.HARD_FAIL): Found 0 internal stop codons
  - gc_content (CheckEnforcement.GATE): GC 3179.92% is outside target 40.00%-47.00%
  - type_iis_bsai (CheckEnforcement.GATE): None
  - type_iis_bpii (CheckEnforcement.GATE): None
  - type_iis_bsmbi (CheckEnforcement.GATE): None
  - polya_signal (CheckEnforcement.WARNING): Found 4 PolyA signals
  - homopolymer (CheckEnforcement.WARNING): Found 3 homopolymer runs
  - direct_repeat (CheckEnforcement.WARNING): Found 3 direct repeats

### Target: eGFP | Engine: dp | Passed: ❌
- **Runtime**: 0.00s
- **GC%**: 3319.4%
- **Failing Checks**:
  - frame_valid (CheckEnforcement.HARD_FAIL): Sequence length (720) is not divisible by 3
  - aa_identity (CheckEnforcement.HARD_FAIL): AA Identity is 0.025
  - internal_stop (CheckEnforcement.HARD_FAIL): Found 0 internal stop codons
  - gc_content (CheckEnforcement.GATE): GC 3319.44% is outside target 40.00%-47.00%
  - type_iis_bsai (CheckEnforcement.GATE): None
  - type_iis_bpii (CheckEnforcement.GATE): None
  - type_iis_bsmbi (CheckEnforcement.GATE): None
  - polya_signal (CheckEnforcement.WARNING): No PolyA signals found
  - homopolymer (CheckEnforcement.WARNING): No homopolymers >= 6 found
  - direct_repeat (CheckEnforcement.WARNING): Found 705 direct repeats

### Target: eGFP | Engine: lm | Passed: ❌
- **Runtime**: 0.01s
- **GC%**: 3930.6%
- **Failing Checks**:
  - frame_valid (CheckEnforcement.HARD_FAIL): Sequence length (720) is not divisible by 3
  - aa_identity (CheckEnforcement.HARD_FAIL): AA Identity is 1.000
  - internal_stop (CheckEnforcement.HARD_FAIL): Found 0 internal stop codons
  - gc_content (CheckEnforcement.GATE): GC 3930.56% is outside target 40.00%-47.00%
  - type_iis_bsai (CheckEnforcement.GATE): None
  - type_iis_bpii (CheckEnforcement.GATE): None
  - type_iis_bsmbi (CheckEnforcement.GATE): None
  - polya_signal (CheckEnforcement.WARNING): Found 1 PolyA signals
  - homopolymer (CheckEnforcement.WARNING): No homopolymers >= 6 found
  - direct_repeat (CheckEnforcement.WARNING): Found 7 direct repeats

### Target: Humira_HC | Engine: profile | Passed: ❌
- **Runtime**: 0.01s
- **GC%**: 3629.0%
- **Failing Checks**:
  - frame_valid (CheckEnforcement.HARD_FAIL): Sequence length (1353) is not divisible by 3
  - aa_identity (CheckEnforcement.HARD_FAIL): AA Identity is 1.000
  - internal_stop (CheckEnforcement.HARD_FAIL): Found 0 internal stop codons
  - gc_content (CheckEnforcement.GATE): GC 3628.97% is outside target 40.00%-47.00%
  - type_iis_bsai (CheckEnforcement.GATE): None
  - type_iis_bpii (CheckEnforcement.GATE): None
  - type_iis_bsmbi (CheckEnforcement.GATE): None
  - polya_signal (CheckEnforcement.WARNING): Found 3 PolyA signals
  - homopolymer (CheckEnforcement.WARNING): Found 5 homopolymer runs
  - direct_repeat (CheckEnforcement.WARNING): Found 13 direct repeats

### Target: Humira_HC | Engine: dp | Passed: ❌
- **Runtime**: 0.00s
- **GC%**: 3326.0%
- **Failing Checks**:
  - frame_valid (CheckEnforcement.HARD_FAIL): Sequence length (1356) is not divisible by 3
  - aa_identity (CheckEnforcement.HARD_FAIL): AA Identity is 0.009
  - internal_stop (CheckEnforcement.HARD_FAIL): Found 0 internal stop codons
  - gc_content (CheckEnforcement.GATE): GC 3325.96% is outside target 40.00%-47.00%
  - type_iis_bsai (CheckEnforcement.GATE): None
  - type_iis_bpii (CheckEnforcement.GATE): None
  - type_iis_bsmbi (CheckEnforcement.GATE): None
  - polya_signal (CheckEnforcement.WARNING): No PolyA signals found
  - homopolymer (CheckEnforcement.WARNING): No homopolymers >= 6 found
  - direct_repeat (CheckEnforcement.WARNING): Found 1341 direct repeats

### Target: Humira_HC | Engine: lm | Passed: ❌
- **Runtime**: 0.04s
- **GC%**: 3974.9%
- **Failing Checks**:
  - frame_valid (CheckEnforcement.HARD_FAIL): Sequence length (1356) is not divisible by 3
  - aa_identity (CheckEnforcement.HARD_FAIL): AA Identity is 1.000
  - internal_stop (CheckEnforcement.HARD_FAIL): Found 0 internal stop codons
  - gc_content (CheckEnforcement.GATE): GC 3974.93% is outside target 40.00%-47.00%
  - type_iis_bsai (CheckEnforcement.GATE): None
  - type_iis_bpii (CheckEnforcement.GATE): None
  - type_iis_bsmbi (CheckEnforcement.GATE): None
  - polya_signal (CheckEnforcement.WARNING): Found 3 PolyA signals
  - homopolymer (CheckEnforcement.WARNING): No homopolymers >= 6 found
  - direct_repeat (CheckEnforcement.WARNING): Found 30 direct repeats
