# DP v2.2.0 research-candidate card

- ID: `dp_v2_2`; version: `2.2.0`; local Python opt-in only.
- Reuses `DPV211Optimizer` transitions, initiation relaxation, STOP/flank
  scanning and LexRank. Adds a penalty at reviewed one-based pause positions.
- Zero penalty or empty positions preserve base sequences and scores on valid
  inputs with identical configuration; engine identity/annotation metadata differ.
- Weights use `codon_weights.get(c, 1e-4)`, not a universal 1e-4 floor. The
  logarithm alone uses a 1e-6 floor. Ramp preference is a proxy, not elongation rate.
- The stem-cell adapter now uses this core with the existing generic human
  reference and provisional policy. No cell-specific calibration is established.
- RNA observations at 37 C remain post-generation model predictions. Passing
  provisional thresholds does not establish PKR avoidance, folding or expression.
- Initiation GC relaxation is disclosed; configured sequence motifs remain vetoes.
- Existing plant `StructuralPauseDPOptimizer`, web and NB-Fold routes retain their
  prior objective. This addition does not silently migrate those routes.

```python
from factorforge.engines.dp_v2_2 import DPV220Optimizer

engine = DPV220Optimizer(
    pause_positions=[2], structure_reference="synthetic:reviewed",
    pause_penalty_weight=2.0, forbidden_motifs=[],
)
result = engine.optimize("KK", {"AAA": 1.0, "AAG": 0.45}, 0, 1)
```

The registry ID returns `DPV220EngineAdapter` and requires explicit
`codon_weights`. Constructor options go in `engine_settings`; no named host
profile is inferred. `python scripts/verify-dp-v2-2-0.py` runs synthetic engineering
checks. These are not a biological benchmark or full manuscript parity proof.
