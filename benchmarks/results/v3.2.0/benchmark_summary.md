# Benchmark Summary (dataset=formal, mode=formal)

## Provenance (reproducibility)
- FactorForge version: 3.2.0
- registry sha256: 8838408157e3e2913c2eb2de2168b369d89cd8ddff9a8c85c4bbced29532cc06
- benchmark_spec sha256: 70adc5fcd681ad3451bc6bd8f685932448e137ec65d49d304b92807e22dca3fe
- input protein fasta sha256: 35dc8ec1f2adaea2be02ad843dd18ca62053595328fe35e38c04bfa596f188e2
- dataset cds fasta sha256: 714a7155c50fff4240c196389ade1860550eb5b949f3cbbd44406ddd6d6cdb53
- run date: 2026-06-13
- codon_table_id: nbenthamiana_legacy_kazusa_sgn_v101
- codon_table_sha256: ddbd0a41da88109a709bca0304581e29bd0a756e4db1c51809d5002e9b2d5e8c
- codon_table_source_status: legacy_metadata_only
- codon_table_build_path_status: incomplete

> **Codon table note:** The current N. benthamiana codon table is a legacy FactorForge reference labeled as derived from Kazusa CodonUsage Database and SGN genome v1.0.1-era resources. The original authoritative build path for the current JSON codon table is incomplete/not verified. The formal benchmark dataset uses SGN QLD183 v103 records; therefore CAI and codon-usage metrics should be interpreted as scores against the configured FactorForge codon reference, not as a de novo SGN QLD183 v103 codon-usage reconstruction.

> This benchmark evaluates in-silico CDS design quality only. It does not demonstrate improved protein expression, yield, folding, or wet-lab performance.

| method | multi_constraint_pass_rate | biological_pass_rate | assembly_pass_rate | gc_in_range_rate | mean_cai |
|---|---|---|---|---|---|
| factorforge_assembly_friendly | 0.63 | 1.00 | 0.66 | 0.97 | 0.9112 |
| factorforge_balanced | 0.27 | 1.00 | 0.28 | 0.98 | 0.9388 |
| factorforge_gc_target | 0.26 | 1.00 | 0.26 | 1.00 | 0.9506 |
| factorforge_high_cai | 0.00 | 1.00 | 0.42 | 0.00 | 0.7343 |
| greedy_cai | 0.22 | 1.00 | 0.31 | 0.73 | 1.0000 |
| native_reference | 0.00 | 1.00 | 0.37 | 0.00 | 0.7249 |
| random_synonymous | 0.01 | 1.00 | 0.26 | 0.01 | 0.7279 |

_native_reference is a biological reference anchor, not an optimizer._
_greedy_cai is a CAI-focused baseline that does not explicitly optimize GC or assembly constraints._

> **Note on factorforge_high_cai:** This profile optimizes against a golden-set codon reference (high-expression gene subset). General CAI scores are computed against the full N. benthamiana CDS-derived codon table and are not expected to be maximized by this profile. Lower general CAI relative to greedy_cai is expected and does not indicate poor performance.

FactorForge preserves amino-acid identity and avoids invalid CDS outputs while improving multi-constraint in-silico CDS design quality relative to simple synonymous-codon baselines.