# True 4-way Sequence Comparison

## Inputs

- comparison_csv: `C:\Work\eijex\factorforge\experiments\results\alpha_run2\alpha_run2_comparison.csv`
- eval_jsonl: `C:\Work\eijex\factorforge\data\training\run4b_pseudolabels_eval.jsonl`
- output_csv: `C:\Work\eijex\factorforge\experiments\results\alpha_run2\true_4way_comparison.csv`
- proteins_evaluated: 3876

## Source Column Check

- present_columns: protein_id, v2_cai, v2_gc_global, v2_local_gc_max, v2_local_gc_min, v3_amino_acid_identity, v3_cai, v3_cai_delta_pct, v3_forbidden_motif_count, v3_gc_global, v3_homopolymer_count, v3_internal_stop_count, v3_invalid_codon_count, v3_local_gc_max, v3_local_gc_min, v3_repeat_count, v3_validator_pass
- missing_v3_sequence_columns: v3_dna_sequence, v3_sequence, v3_decoded_sequence, v3_decoded_cds, v3_cds
- v3_sequence_status: redecode_required_no_local_sequence

Candidate C uses stored v3 decoded CDS sequences only if such a column exists. The current local alpha_run2 comparison stores v3 metrics but not decoded CDS, so v3 rows are metric placeholders and require re-decode to support codon identity.

## 4-way Comparison Summary

| candidate_type | n | sequence_available | CAI mean | GC mean | repeat mean | homopolymer mean | motif mean | score mean | validator pass |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| v2_high_cai | 3876 | 1.0000 | 0.7340 | 33.2250 | 5.8599 | 2.1734 | 0.0000 | 0.3329 | 1.0000 |
| v2_gc_target | 3876 | 1.0000 | 0.8025 | 42.5441 | 3.7763 | 0.3096 | 0.0000 | 0.6982 | 1.0000 |
| v3_run4b_decoded | 3876 | 0.0000 | 0.8025 | 42.5440 | 3.7763 | 0.3088 | 0.0000 | 0.6983 | 1.0000 |
| feasibility_best | 3876 | 1.0000 | 0.9605 | 54.9231 | 5.2771 | 0.0044 | 0.0000 | 0.7965 | 1.0000 |

## v3 vs v2 gc_target Codon Identity

- mean codon identity: N/A
- interpretation: Cannot determine from local artifacts because v3 decoded CDS sequences are not stored in alpha_run2 outputs. Re-decode with the Run4B checkpoint is required.

## Feasibility Candidate

- feasibility_source_counts: {'dp_optimal_gc_40_55': 3876}
- `dp_optimal_gc_40_55` means `analyze_feasibility()` selected the maximum-CAI sequence under global GC 40-55%.

## Candidate Status Counts

```json
{
  "redecode_required_no_local_sequence": 3876,
  "sequence_available": 11628
}
```

## Honest Conclusion

Cannot determine from local artifacts because v3 decoded CDS sequences are not stored in alpha_run2 outputs. Re-decode with the Run4B checkpoint is required.

This is an in-silico sequence comparison. It does not change v2 production behavior and does not support yield or expression claims.
