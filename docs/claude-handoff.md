# Claude Code Handoff Template

Use this template when a CI failure or development idea needs a focused Claude
Code review.

## Request

```text
Repository:
Issue:
Branch:
Commit:
Workflow:
Run URL:

Task:
Review the failure or idea and propose the smallest safe change.

Scope:
- Inspect the relevant workflow, source, tests, and docs.
- Do not auto-merge or push changes without maintainer approval.
- Preserve the in-silico validation language in public documentation.
- Do not open public security vulnerability details.

Done when:
- Root cause or implementation path is documented.
- Patch is prepared if appropriate.
- Proportional verification is run or explicitly marked N/A.
```

## Notes

- For CI failures, include the failed job name and the first relevant log lines.
- For Docker failures, include the image target and registry step.
- For docs failures, include the page path and MkDocs output.
- For wet-lab feedback, summarize only non-sensitive metadata unless explicit
  permission exists to publish details.
