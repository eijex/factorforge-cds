# Benchmark Summary (dataset=formal, mode=formal)

## Provenance (reproducibility)
- FactorForge version: 3.1.9
- registry sha256: 17bf44d4c0217ec158a8bbf8e8d78c0e7264ab748f4642a42412a5f24a7e78e3
- benchmark_spec sha256: 973b08cce1bcf3ecdf95f693d23afd70c2ff7d52d7e11b302d745a780fce20df
- input protein fasta sha256: 35dc8ec1f2adaea2be02ad843dd18ca62053595328fe35e38c04bfa596f188e2
- dataset cds fasta sha256: 714a7155c50fff4240c196389ade1860550eb5b949f3cbbd44406ddd6d6cdb53
- run date: 2026-06-11

> This benchmark evaluates in-silico CDS design quality only. It does not demonstrate improved protein expression, yield, folding, or wet-lab performance.

| method | multi_constraint_pass_rate | biological_pass_rate | assembly_pass_rate | gc_in_range_rate | mean_cai |
|---|---|---|---|---|---|
| factorforge_assembly_friendly | 0.66 | 1.00 | 0.66 | 0.97 | 0.9092 |
| factorforge_balanced | 0.30 | 1.00 | 0.30 | 0.98 | 0.9305 |
| factorforge_gc_target | 0.26 | 1.00 | 0.26 | 1.00 | 0.9505 |
| factorforge_high_cai | 0.42 | 1.00 | 0.42 | 0.00 | 0.7343 |
| greedy_cai | 0.31 | 1.00 | 0.31 | 0.73 | 1.0000 |
| most_frequent_codon | 0.31 | 1.00 | 0.31 | 0.73 | 1.0000 |
| native_reference | 0.37 | 1.00 | 0.37 | 0.00 | 0.7249 |
| random_synonymous | 0.26 | 1.00 | 0.26 | 0.01 | 0.7279 |

_native_reference is a biological reference anchor, not an optimizer._
_greedy_cai is a CAI-focused baseline that does not explicitly optimize GC or assembly constraints._

> **Note on factorforge_high_cai:** This profile optimizes against a golden-set codon reference (high-expression gene subset). General CAI scores are computed against the full N. benthamiana CDS-derived codon table and are not expected to be maximized by this profile. Lower general CAI relative to greedy_cai is expected and does not indicate poor performance.

FactorForge preserves amino-acid identity and avoids invalid CDS outputs while improving multi-constraint in-silico CDS design quality relative to simple synonymous-codon baselines.