# Inbox Setup

The failure inbox captures failed GitHub Actions runs as GitHub issues and sends
an optional Discord notification.

## GitHub Secret

Create a repository secret named `DISCORD_WEBHOOK_URL`:

1. Open repository settings.
2. Go to **Secrets and variables** -> **Actions**.
3. Add `DISCORD_WEBHOOK_URL` with the Discord webhook for `#ci-failures`.

If the secret is missing or Discord rejects the request, the workflow skips or
continues without failing the inbox run.

## Labels

The workflow creates these labels if they do not already exist:

- `ci-failure`
- `claude-review`
- `docker`
- `tests`
- `lint`
- `deploy`
- `docs`
- `unknown`
- `urgent`
- `idea`
- `later`
- `needs-triage`
- `github-actions`
- `user-feedback`
- `feature-request`
- `bug-report`
- `docs-feedback`
- `wet-lab`

Deploy failures also receive `urgent`.

## Classification

Workflow names are classified by keyword:

| Keywords | Category |
| --- | --- |
| `docker`, `build`, `container`, `image` | `docker` |
| `test`, `pytest`, `unit` | `tests` |
| `lint`, `ruff`, `format` | `lint` |
| `deploy`, `vercel`, `pages` | `deploy` |
| `docs`, `mkdocs` | `docs` |
| Any other name | `unknown` |

## Safety Policy

- The inbox creates issues and notifications only; it does not change code.
- Do not auto-merge, auto-push, or auto-open pull requests from inbox issues.
- Do not publish secrets, private collaborator data, or vulnerability details in
  public issues.
- Keep public claims aligned with `VALIDATION.md` and `docs/validation.md`.
