# Codon Table Provenance — N. benthamiana

## Overview

FactorForge uses a pre-built codon usage table for *Nicotiana benthamiana* (`src/factorforge/data/nbenthamiana_codons.json`) to drive codon optimization scoring (CAI) and sequence generation.

This document discloses the provenance status of that table and the limitations users should be aware of when interpreting results.

---

## Current Production Codon Reference

| Field | Value |
|-------|-------|
| **ID** | `nbenthamiana_nbev11_hc_v2` |
| **Species** | *Nicotiana benthamiana* |
| **NCBI TaxID** | 4100 |
| **File** | `src/factorforge/data/nbenthamiana_codons.json` |
| **SHA-256** | See `data/reference/active_codon_reference.json` and `data/reference/codon_table_manifest_nbev11_hc_v2.json` |
| **Source label** | NbeV1.1 high-confidence CDS-derived native-composition anchor |
| **Source status** | `official_packaged` / active production software default |
| **Build path status** | Manifested and checksum-guarded |
| **Authoritative build script** | `scripts/build_codon_profile.py` with `strict_nuclear_cds_v1` filtering |

---

## Provenance Disclosure

The current *N. benthamiana* production software default is `nbenthamiana_nbev11_hc_v2`, an NbeV1.1 high-confidence CDS-derived codon-reference asset with the default GC policy band aligned to 40-47%. The legacy Kazusa/SGN v1.0.1-era table is retained as historical comparator/provenance data, not the public product default. The formal benchmark dataset uses SGN QLD183 v103 records; therefore CAI and codon-usage metrics should be interpreted as scores against the configured FactorForge codon reference, not as comparative biological-performance evidence.

---

## Known Limitations

1. The current production default is an in-silico software reference, not wet-lab validation, expression-output prediction, or comparative biological-performance evidence.
2. The legacy Kazusa/SGN composite remains packaged for historical continuity but is not the public product default.
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

## v3.3.0 Default Activation

FactorForge v3.3.0 activates `nbenthamiana_nbev11_hc_v2` as the current *N. benthamiana* production software default:

- **Active reference**: NbeV1.1 high-confidence CDS-derived codon usage
- **Default GC policy band**: 40-47%
- **Policy**: Public API/UI codon-reference overrides remain unsupported
- **Comparator assets**: Legacy Kazusa/SGN, NbeV1.1 all-CDS, QLD183 v103, and other research/candidate assets remain retained for provenance and controlled analysis only

Current production results should be interpreted against the active `nbenthamiana_nbev11_hc_v2` software reference unless a controlled research/benchmark path explicitly pins a comparator table.

---

## Machine-Readable Manifest

See [`data/reference/codon_table_manifest.json`](../../data/reference/codon_table_manifest.json) for the machine-readable version of this provenance record, validated against [`schemas/codon_table_manifest.schema.json`](../../schemas/codon_table_manifest.schema.json).
