# v3-alpha Evidence/Reranking Evaluation

## Inputs

- comparison_csv: `C:\Work\eijex\factorforge\experiments\results\run4b_full\run4_comparison.csv`
- eval_jsonl: `C:\Work\eijex\factorforge\data\training\run4b_pseudolabels_eval.jsonl`
- extended_csv: `C:\Work\eijex\factorforge\experiments\results\run4b_full\run4_comparison_extended.csv`
- proteins_evaluated: 3876

## Source CSV Columns

- present: protein_id, v2_cai, v2_gc_global, v2_local_gc_max, v2_local_gc_min, v3_amino_acid_identity, v3_cai, v3_cai_delta_pct, v3_forbidden_motif_count, v3_gc_global, v3_homopolymer_count, v3_internal_stop_count, v3_invalid_codon_count, v3_local_gc_max, v3_local_gc_min, v3_repeat_count, v3_validator_pass
- missing: v2_forbidden_motif_count, v2_homopolymer_count, v2_repeat_count, v2_sequence, v3_sequence

The Run 4B comparison CSV does not include decoded CDS sequence columns. The v3 Run4B row therefore uses the decoded metrics already recorded in `run4_comparison.csv`; v2 high_cai, v2 gc_target, and feasibility_best rows are recomputed from sequence-level APIs.

## 4-way Comparison Summary

| candidate_type | n | score_mean | score_median | score_min | score_max | cai_mean | gc_mean | gc_in_range_rate | validator_pass_rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| feasibility_best | 3876 | 0.7965 | 0.8202 | -1.4835 | 1.1000 | 0.9605 | 54.9231 | 1.0000 | 1.0000 |
| v2_gc_target | 3876 | 0.6982 | 0.7234 | -1.3171 | 0.9903 | 0.8025 | 42.5441 | 1.0000 | 1.0000 |
| v2_high_cai | 3876 | 0.3329 | 0.3833 | -2.6911 | 0.8007 | 0.7340 | 33.2250 | 0.0052 | 1.0000 |
| v3_run4b_decoded | 3876 | 0.6983 | 0.7235 | -1.3171 | 0.9903 | 0.8025 | 42.5440 | 1.0000 | 1.0000 |

## Multi-objective Score

Score formula: `cai + 0.1 * GC_in_40_55 - 0.05 * repeats - 0.05 * homopolymers - 0.1 * forbidden_motifs`.

## v3 Advantage Cases

- v3 score > v2 high_cai: 3843 / 3876
- v3 score > v2 gc_target: 1931 / 3876
- first v3 > v2 high_cai IDs: ['NbL04g19790.1', 'NbL01g19020.1', 'NbL12g19450.1', 'NbL13g12940.1', 'NbL06g16690.1', 'NbL18g13520.1', 'NbL16g17570.1', 'NbL14g12980.1', 'NbL03g24230.1', 'NbL05g22170.1', 'NbL02g16680.1', 'NbL05g01090.1', 'NbL04g06350.1', 'NbL06g05440.1', 'NbL14g17400.1', 'NbL16g14100.1', 'NbL13g17530.1', 'NbL10g01790.1', 'NbL15g22530.1', 'NbL14g00650.1']
- first v3 > v2 gc_target IDs: ['NbL18g13520.1', 'NbL14g12980.1', 'NbL03g24230.1', 'NbL02g16680.1', 'NbL05g01090.1', 'NbL04g06350.1', 'NbL16g14100.1', 'NbL10g01790.1', 'NbL19g08970.1', 'NbL04g06940.1', 'NbL06g17910.1', 'NbL13g07490.1', 'NbL03g16170.1', 'NbL11g20730.1', 'NbL05g05610.1', 'NbL03g17710.1', 'NbL14g20930.1', 'NbL11g08560.1', 'NbL12g20570.1', 'NbL10g16230.1']

## Protein Characteristics Associated With v3 Advantage

v3 > v2 high_cai:

```json
{
  "count": 3843,
  "protein_length_mean": 329.6247723132969,
  "selected_teacher_type_counts": {
    "gc_target": 3841,
    "high_cai": 2
  }
}
```

v3 > v2 gc_target:

```json
{
  "count": 1931,
  "protein_length_mean": 328.3542206110823,
  "selected_teacher_type_counts": {
    "gc_target": 1928,
    "high_cai": 3
  }
}
```

## Honest Conclusion

v3 is comparable

This is an in-silico candidate evidence comparison only. It does not change v2 production behavior and does not support yield or expression claims.
