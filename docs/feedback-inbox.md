# Feedback Inbox

This document describes the planned v0.2 feedback path. It is documentation
only; no Google Form, Make scenario, or GitHub automation is implemented here.

## Pipeline

```text
Google Form -> Make -> GitHub Issue -> Label triage -> Maintainer review
```

The form is the preferred entry point for user feedback because it can collect
structured fields before creating an issue.

## Issue Labels

- `user-feedback`: general feedback entry
- `bug-report`: reported unexpected behavior
- `feature-request`: requested capability
- `docs-feedback`: documentation issue or suggestion
- `wet-lab`: experimental feedback requiring careful review
- `joss-review`: publication review follow-up
- `needs-triage`: maintainer triage needed

## Don Stewart Wet-Lab Data Schema

Use this schema for structured experimental feedback. Keep private or unpublished
details out of public issues unless the submitter has explicitly approved
publication.

| Field | Purpose |
| --- | --- |
| `submitter_name` | Person submitting the result |
| `submitter_contact` | Email or preferred contact path |
| `construct_id` | Internal or public construct identifier |
| `host_system` | Host organism or expression system |
| `target_protein` | Protein or construct target |
| `factorforge_version` | Version used to generate the CDS |
| `input_sequence_hash` | Hash or non-sensitive reference to input sequence |
| `output_sequence_hash` | Hash or non-sensitive reference to output sequence |
| `assay_type` | Assay or measurement type |
| `observed_result` | Structured observation, not a broad performance claim |
| `controls` | Control construct or baseline notes |
| `replicate_count` | Number of replicates |
| `raw_data_location` | Private storage pointer when data cannot be public |
| `publication_permission` | Whether public issue details are allowed |
| `notes` | Additional non-sensitive context |

## Google Opal UI

Google Opal can be used as an optional front-end for guided feedback capture:

- Step users through category selection.
- Collect only the fields needed for the selected category.
- Show a privacy reminder before submission.
- Send the final structured payload to the Make scenario.

GitHub issues remain the tracking system after submission.
