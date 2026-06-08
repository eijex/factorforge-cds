# Feedback Inbox

This document describes a planned feedback path. It is documentation only; no Google Form, Make scenario, or GitHub automation is implemented here.

## Pipeline

```text
Google Form -> maintainer review -> public-safe GitHub Issue or documentation update
```

The form may collect structured public-safe summary fields before maintainer review. Public GitHub Issues must not contain private or sensitive wet-lab data.

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

## Google Opal UI

Google Opal can be used as an optional front-end for guided feedback capture:

- Step users through category selection.
- Collect only public-safe fields needed for the selected category.
- Show a privacy reminder before submission.
- Send the final structured payload to maintainer review.

GitHub issues remain public and should contain only reviewed public-safe summaries.
