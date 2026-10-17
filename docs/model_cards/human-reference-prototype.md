# Experimental human-reference prototype 1.0.1

This local interface retains historical StemCell names for import compatibility. Human-frequency table provenance is attributed to a legacy Kazusa source, with cell-specific calibration unverified; values are unchanged. SCI criteria and approval have not been obtained.

The DP objective remains codon log-adaptiveness minus an optional reviewed-position penalty, with GC and motif constraints. With no positions supplied this is sequence-only, not a protein-folding-aware design. Explicit positions require positions_reviewed=True and a structure_reference. Neither label is a server authorization or biological validation.

ViennaRNA uses isolated 37 C model details. Full sense-CDS, first 60 nt and complete 40-nt windows at 5-nt steps are modeled after generation. No UTR/vector/actual transcript is inferred. MFE is a modeled observation, not measured initiation, folding kinetics or expression. No 5-prime relaxation or RNA-constrained DP pruning is performed.

The packaged policy is mirrored at configs/sops/human_stemcell_conservative_reference.json and validated against known rules. Legacy initiation rule is ignored because its -3 threshold/window is a different contract; prototype observation checks explicitly use the existing PRD's provisional first-60-nt minimum -5 and 40-nt minimum -12 kcal/mol and stem-limit comparison <20. These numbers are not calibrated efficacy/safety thresholds. Existing GC defaults 44–50 percent and codon weights are preserved. miRNA-related strings are exact-pattern exclusions only, not a complete target predictor or proof of silencing avoidance.

Computational status is PASS, BLOCKED or INCOMPLETE. Overall status remains REVIEW_REQUIRED even if computational checks pass; biological_validation and pkr_assessment remain NOT_AVAILABLE. Missing ViennaRNA or invalid output cannot produce a computational pass. A failing post-filter retains the original candidate for review and does not silently regenerate or claim a globally RNA-optimal solution.

Example on a synthetic protein:

```python
from factorforge.engines.stemcell_3dp import StemCell3DPOptimizer
result = StemCell3DPOptimizer(target_gc_min=0, target_gc_max=1).optimize("MASTWQHLEK")
assert result["status"] == "REVIEW_REQUIRED"
# Bounds above are a synthetic demonstration, not a laboratory policy.
```

Evaluate a supplied sense CDS with evaluate_candidate(dna, expected_protein) to use the same checks without modifying it. Preserve source, exact protein/construct identity and configuration alongside any candidate comparison. Human MSC and iPSC experiments require separate reviewed contexts. No new host/API/deployment, yield model or approved experimental protocol is provided.
