# ValidationHub — Prohibited Fields

This document defines fields that must not appear in any ValidationHub intake record.

## Policy

ValidationHub intake records are designed to collect only non-confidential, non-identifiable information. The following field names are explicitly prohibited to protect:

1. Confidentiality of proprietary or unpublished sequences
2. Privacy of submitters and collaborators
3. Integrity of the public validation record

## Prohibited Field Names

### Sequence Fields

| Field name | Reason |
|-----------|--------|
| `raw_sequence` | Direct sequence disclosure |
| `sequence` | Direct sequence disclosure |
| `dna_sequence` | Direct sequence disclosure |
| `rna_sequence` | Direct sequence disclosure |
| `protein_sequence` | Direct sequence disclosure |
| `amino_acid_sequence` | Direct sequence disclosure |
| `construct_sequence` | Direct sequence disclosure |
| `customer_sequence` | Direct sequence disclosure |
| `proprietary_sequence` | Direct sequence disclosure |
| `sequence_hash` | Indirect sequence proxy — can enable de-anonymization of private sequences in intake context |

### Identifier Fields

| Field name | Reason |
|-----------|--------|
| `construct_id` | Internal identifier, potentially linkable to proprietary data |
| `customer_construct_id` | Links intake to customer-specific assets |
| `internal_batch_id` | Internal production identifier |
| `protein_name` | May reveal proprietary target |

### Contact / Personal Fields

| Field name | Reason |
|-----------|--------|
| `submitter_name` | PII |
| `submitter_email` | PII |
| `respondent_email` | PII |
| `contact_email` | PII |
| `phone_number` | PII |
| `personal_identifier` | PII |

## Clarification: sequence_hash in Design Package

`sequence_hash` is listed above as prohibited in **intake records**.

It is allowed in the **Design Package** (`design_package.json` → `evidence.sequence_hash`), because the Design Package is a self-contained computational artifact — not an intake record submitted by a user. These are different contexts.

## Using Public Sequences

For publicly citable sequences (e.g., sfGFP from PDB:2B3P), do **not** include the raw sequence in the intake record. Instead, use:

```json
"sequence_publication_source": "Pédelacq et al. 2006, Nat Biotechnol 24:79-88, PMID 16369541"
```

The `design_package_id` field links the intake record back to the corresponding Design Package, which contains the `evidence.sequence_hash` for verification purposes.

## Prohibited Submission Types

Do not submit intake records based on:

- Private, proprietary, or unpublished sequences
- Patient-derived sequences
- Customer or collaborator sequences without explicit consent
- Internal batch or development sequences

> "Do not submit private, proprietary, unpublished, patient-derived, customer, collaborator, or internal raw sequences through public ValidationHub records. Use `sequence_publication_source` for publicly citable sequences."

## Source of Truth

The canonical prohibited field list is maintained in the ValidationHub private repository:
`validationHub/schemas/prohibited_fields.json`

This document is a human-readable summary derived from that source.
