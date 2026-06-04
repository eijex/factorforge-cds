# Developer Inbox Roadmap

The developer inbox is a planned workflow for turning operational signals and
user feedback into reviewable implementation tasks. It is not an auto-fix
system.

## v0.1: GitHub Actions Failure Inbox

Status: planned, not implemented.

The first version will collect CI failures and create a structured review path:

- GitHub Actions `workflow_run` event on failure creates a GitHub Issue
  automatically.
- Discord webhook notification alerts maintainers.
- Claude Code handoff document includes relevant files, likely scope, and test
  commands.
- Safety rule: no auto-fix, no auto-PR, no auto-merge, no auto-deploy.

Planned files:

- `.github/workflows/failure-inbox.yml`
- `.github/ISSUE_TEMPLATE/idea.yml`

The failure inbox should package context for a developer or agent-assisted
coding session. It should not mutate code or scoring behavior.

## v0.2: User Feedback Inbox

Status: planned, not implemented.

The second version will collect user feedback and convert it into triage-ready
issues:

- Google Form or email feedback creates a GitHub Issue.
- Discord alert notifies maintainers.
- Claude triage classification labels the item for review.

Planned files:

- `.github/ISSUE_TEMPLATE/feedback.yml`
- `docs/feedback-inbox.md`

Feedback triage should keep bug reports, feature requests, documentation
requests, and biological-design questions separate. Biological or validation
claims should require evidence review before being accepted into product docs.
