# DP v2.1 development-candidate card

| Field | Value |
|---|---|
| Engine ID | `dp_v2_1` |
| Engine version | `2.1.0-dev` |
| Product line | FactorForge `3.5.0` release candidate |
| Status | Development candidate; explicit opt-in |
| Current host | *Nicotiana benthamiana* |
| Default replacement | No; `feasibility_best` continues to use DP v2 `2.0.1` |

DP v2.1 performs deterministic synonymous CDS candidate generation across
three declared axes:

1. `HARD` — requested global GC range and configured Type IIS motif avoidance.
2. `OPTIMIZED` — host-relative codon adaptation.
3. `OPTIMIZED` — position-dependent 5′ initiation scoring using an
   open-topology codon proxy over the initial ramp.

The third axis is not an RNA-folding calculation. MFE is returned as
`not_computed`; any folding result belongs to a separate
`INDEPENDENTLY_EVALUATED` stage.

Outputs are in-silico design candidates. This card does not claim improved
expression, yield, folding, synthesis acceptance, or wet-lab validation.
