# Worked Example: sfGFP CDS Design for *N. benthamiana*

A fully reproducible end-to-end CDS design example using Superfolder GFP (sfGFP) as the target protein.

## Overview

| Parameter | Value |
|-----------|-------|
| Input protein | sfGFP (Pédelacq et al. 2006, 236 aa, PDB:2B3P) |
| Optimization profile | `assembly_friendly` |
| Host | *Nicotiana benthamiana* |
| Seed | 320 (deterministic) |
| Scoring contract | v1.1 |

## Prerequisites

```bash
pip install factorforge-cds
```

## Input: sfGFP Protein Sequence

The input is the Superfolder GFP amino acid sequence from Pédelacq et al. 2006, obtained from PDB entry 2B3P (chain A). The chromophore residue (X in the PDB ATOM record) is substituted with Tyr (Y), as this is the pre-cyclization amino acid. The C-terminal His-tag is removed.

**Citation**: Pédelacq, J.-D., Cabantous, S., Tran, T., Terwilliger, T. C. & Waldo, G. S. Engineering and characterization of a superfolder green fluorescent protein. *Nat Biotechnol* **24**, 79–88 (2006). PMID: 16369541

The sequence is stored in `input_sequence.faa`.

## Run

```bash
# From this directory:
python run_example.py
```

Expected output (if frozen outputs already exist):
```
OK: frozen outputs match current run (reproducible)
```

To regenerate frozen outputs from scratch:
```bash
python run_example.py --freeze
```

Expected output:
```
Frozen: .../output/design_package.json
Frozen: .../output/validation_summary.json

Key results:
  CAI:                 0.769
  GC%:                 56.78
  gc_in_target_range:  True
  multi_constraint_pass: True
```

## Output Files

### `output/design_package.json`

A FactorForge Open Bio Design Package (schema v1.0). Key fields:

| Field | Value | Description |
|-------|-------|-------------|
| `metrics.cai` | 0.769 | Codon Adaptation Index |
| `metrics.gc_percent` | 56.78 | GC content (%) |
| `validation.biological_pass` | true | AA identity = 1.0, no internal stops |
| `validation.assembly_pass` | true | No forbidden Golden Gate TypeIIS sites |
| `claim_boundary.in_silico_only` | true | Computational prediction only |

### `output/validation_summary.json`

A ValidationHub computational intake record (schema v0.1). Key fields:

| Field | Value | Description |
|-------|-------|-------------|
| `evidence_type` | `computational_validation` | No wet-lab data |
| `computational.scoring_contract_version` | `v1.1` | Contract version |
| `computational.gc_in_target_range` | true | 55–65% target range |
| `computational.multi_constraint_pass` | true | Derived: bio AND assembly AND gc_in_range |
| `experimental.outcome_status` | `not_tested` | Wet-lab not yet performed |

## Reproducibility Note

**Scoring contract v1.1** defines `multi_constraint_pass` as:

```
biological_pass AND assembly_pass AND gc_in_target_range
```

`multi_constraint_pass` is always a **derived** field — it is never asserted directly. The primitives (`biological_pass`, `assembly_pass`, `gc_in_target_range`) are computed independently from the optimized DNA sequence using production FactorForge validators.

## ValidationHub Connection

`validation_summary.json` conforms to the ValidationHub computational intake schema (`docs/validationhub/intake_schema_v0.1.public.json`). The `design_package_id` field links this record back to the corresponding `design_package.json`.

When wet-lab expression data becomes available, the `experimental` section can be completed:

```json
"experimental": {
  "assay_type": "western_blot",
  "outcome_status": "pass"
}
```

## Citation

If you use this worked example in publications, please cite:

```
Pédelacq, J.-D. et al. Engineering and characterization of a superfolder green
fluorescent protein. Nat Biotechnol 24, 79–88 (2006). https://doi.org/10.1038/nbt1172
```
