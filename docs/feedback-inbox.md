# Feedback Inbox

This document describes a planned feedback path. It is documentation only; no external form automation is implemented here.

## Submission Channels

| Channel | Use |
| --- | --- |
| [GitHub Issue (wet_lab_result)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) | Public-safe coarse summaries only |
| [Google Form](https://docs.google.com/forms/d/e/1FAIpQLSeSx-wYvF6YwHhSPdLMl-L44frCugdm25X_eDz50OaqTD66qA/viewform) | Structured wet-lab result submission |
| Email: eijex.lab@gmail.com | Private or sensitive summaries |

Public GitHub Issues must not contain private or sensitive wet-lab data.

## Issue Labels

- `user-feedback`: general feedback entry
- `bug-report`: reported unexpected behavior
- `feature-request`: requested capability
- `docs-feedback`: documentation issue or suggestion
- `wet-lab`: public-safe experimental feedback requiring careful review
- `needs-triage`: maintainer triage needed

## Public-Safe Wet-Lab Summary Fields

Use this shape for public experimental feedback. Keep private, unpublished, identifying, or confidential details out of public issues.

| Field | Purpose |
| --- | --- |
| `factorforge_version` | Version used to generate the CDS |
| `optimization_profile` | Public optimization profile name |
| `host_organism` | Host organism or expression system |
| `protein_class` | Coarse non-confidential protein class |
| `expression_system_assay_type` | Assay or measurement type |
| `harvest_timepoint` | Coarse harvest timepoint |
| `replicate_category` | Coarse replicate category |
| `comparison_result` | Improved, equivalent, reduced, inconclusive, no control, or not reported |
| `expression_result` | Detected, not detected, weak, strong, inconclusive, or not reported |
| `approximate_yield_range` | Coarse yield range, not exact confidential process data |
| `institution_disclosure` | Anonymous, disclosed, or not provided |
| `non_confidential_notes` | Additional non-sensitive context |

Do not include submitter identity, private contact information, raw sequences, sequence hashes, internal batch IDs, confidential construct identifiers, private customer names, patient data, exact confidential process parameters, or unpublished construct-identifying labels.

GitHub issues remain public and should contain only public-safe summaries.
