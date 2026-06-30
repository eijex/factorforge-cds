# Codon Table Provenance — N. benthamiana

## Overview

FactorForge uses a pre-built codon usage table for *Nicotiana benthamiana* (`src/factorforge/data/nbenthamiana_codons.json`) to drive codon optimization scoring (CAI) and sequence generation.

This document discloses the provenance status of that table and the limitations users should be aware of when interpreting results.

---

## Current Codon Table

| Field | Value |
|-------|-------|
| **ID** | `nbenthamiana_legacy_kazusa_sgn_v101` |
| **Species** | *Nicotiana benthamiana* |
| **NCBI TaxID** | 4100 |
| **File** | `src/factorforge/data/nbenthamiana_codons.json` |
| **SHA-256** | `ddbd0a41da88109a709bca0304581e29bd0a756e4db1c51809d5002e9b2d5e8c` |
| **Source label** | Kazusa CodonUsage Database + SGN genome v1.0.1 |
| **Source status** | `legacy_metadata_only` |
| **Build path status** | `incomplete` |
| **Authoritative build script** | Not available |

---

## Provenance Disclosure

The current N. benthamiana codon table is a legacy FactorForge reference labeled as derived from Kazusa CodonUsage Database and SGN genome v1.0.1-era resources. The original authoritative build path for the current JSON codon table is incomplete/not verified. The formal benchmark dataset uses SGN QLD183 v103 records; therefore CAI and codon-usage metrics should be interpreted as scores against the configured FactorForge codon reference, not as a de novo SGN QLD183 v103 codon-usage reconstruction.

---

## Known Limitations

1. The original script or exact command used to generate the current JSON codon table is not available or not verified.
2. The current codon table is not a freshly rebuilt SGN QLD183 v103 codon table.
3. Benchmark CAI scores should be interpreted against the configured FactorForge codon reference.

---

## Archive Script Disclaimer

The repository contains a historical script at:

```
archive/v3-ml-prototype/scripts/1_data_preparation/build_codon_table.py
```

**This script is NOT the authoritative builder for the current production codon table.**

It was developed as part of the v3 ML prototype and requires expression-level input (Salmon `quant.sf` / TPM values) to weight codon frequencies. Its output format (CSV) also differs from the current production JSON format. It is preserved for historical reference only and must not be cited as the provenance source for `nbenthamiana_codons.json`.

---

## Benchmark Relationship

The v3.2.0 formal benchmark uses:

- **Benchmark dataset**: SGN QLD183 v103 CDS records (N=49,257)
- **Codon table**: `nbenthamiana_nbev11_hc_v2` (NbeV1.1 high-confidence CDS-derived current software default); legacy `nbenthamiana_legacy_kazusa_sgn_v101` retained as historical comparator

Using a codon table derived from one genome version while benchmarking on another is standard practice in codon optimization. The codon table is the *scoring reference*, not the benchmark dataset itself. CAI values reflect how well each sequence matches the configured FactorForge reference.

---

## Claim Boundary

This codon table is a configured in-silico FactorForge reference for codon-usage scoring and optimization. It is not wet-lab validation, yield prediction, or expression-success evidence.

---

## Planned Rebuild — v3.3.0

A full rebuild of the N. benthamiana codon table from SGN QLD183 v103 is planned for FactorForge v3.3.0:

- **Target reference**: SGN QLD183 v103
- **Approach**: Download CDS with documented script, apply explicit filtering rules, compute codon frequencies
- **Expected effect**: Codon weights, CAI scores, profile rankings, and benchmark results may change
- **Deliverables**: New `codon_table_id`, updated manifest, old vs. new codon weight diff, re-run benchmark (seed=320)

Current production results should be interpreted against the active `nbenthamiana_nbev11_hc_v2` software reference unless a controlled research/benchmark path explicitly pins a comparator table.

---

## Machine-Readable Manifest

See [`data/reference/codon_table_manifest.json`](../../data/reference/codon_table_manifest.json) for the machine-readable version of this provenance record, validated against [`schemas/codon_table_manifest.schema.json`](../../schemas/codon_table_manifest.schema.json).
