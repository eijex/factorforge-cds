# Feedback Inbox

This document describes a planned feedback path. It is documentation only; no external form automation is implemented here.

## Submission Channels

| Channel | Use |
| --- | --- |
| [Share Wet-lab Results (GitHub)](https://github.com/eijex/factorforge-cds/issues/new?template=wet_lab_result.yml) | Public-safe coarse summaries only |
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
| `version` | FactorForge version used to generate the CDS |
| `profile` | Public optimization profile used |
| `host` | Host context |
| `protein_class` | Coarse non-confidential protein class |
| `assay_category` | Coarse assay category |
| `comparison_result` | Improved, equivalent, reduced, inconclusive, no control, or not reported |
| `expression_result` | Signal observed, no signal observed, inconclusive, not measured, or not reported |
| `replicate_category` | Coarse replicate category |
| `non_confidential_notes` | Additional non-sensitive context |

Do not include submitter identity, private contact information, raw sequences, sequence hashes, internal batch IDs, confidential construct identifiers, private customer names, patient data, exact confidential process parameters, or unpublished construct-identifying labels.

GitHub issues remain public and should contain only public-safe summaries.
