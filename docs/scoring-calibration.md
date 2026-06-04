# Scoring Calibration

FactorForge separates optimization targets, feasibility lower bounds, and repair
guardrails. These values should not be treated as interchangeable.

## CAI Thresholds

- `target_cai`: The optimizer aims for this value during candidate generation.
  A typical target is `0.90`.
- `min_feasible_cai`: The dynamic-programming feasibility engine lower bound.
  A typical lower bound is `0.82`.
- `repair_floor_cai`: The dinucleotide repair guard floor. A typical guard
  floor is `0.75`.

`target_cai` describes the preferred design goal. `min_feasible_cai` describes
whether a candidate space can reasonably satisfy a requested design. The repair
floor protects post-generation fixes from preserving a motif at the cost of
excessive codon-adaptation loss.

## GC Band Scoring

The old point-distance formula was:

```python
1.0 - abs(gc - gc_opt) / 50.0
```

That formula was too weak outside the desired target because distant GC values
could retain a relatively high score.

The current band function is:

```python
def gc_band_score(
    gc: float,
    gc_min: float,
    gc_max: float,
    decay_width: float = 20.0,
) -> float:
```

The function returns a full score inside `[gc_min, gc_max]`. Outside the band,
the score decays linearly to `0.0` over `decay_width` percentage points. For the
current plant default, the accepted band is `55.0` to `65.0` GC with a
`20.0`-point decay width.

## Composite Score Weights

The live scoring configuration currently defines these profile weights:

| Profile | w_cai | w_gc | w_mfe | Notes |
|---------|------:|-----:|------:|-------|
| `balanced` | 0.5 | 0.3 | 0.2 | Default balanced optimization |
| `high_cai` | 0.8 | 0.1 | 0.1 | Strong CAI pressure |
| `gc_target` | 0.1 | 0.7 | 0.2 | GC target priority |
| `assembly_friendly` | 0.3 | 0.4 | 0.3 | Reduced CAI pressure, higher GC/MFE weight |

Weights are normalized by `ScoringConfig` over active components. MFE scoring is
used only when sequence context is available and ViennaRNA bindings are
available. Otherwise MFE contribution is disabled for that calculation.

Future host-specific scoring configuration should live in host profile YAML.
That migration should include snapshot tests so changes to CAI, GC, MFE, and
dinucleotide semantics are visible before release.
