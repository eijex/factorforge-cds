# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 3.1.x (latest) | ✅ Security patches applied |
| 3.0.x | ❌ No longer supported |
| < 3.0 | ❌ No longer supported |

Only the latest patch release of the current minor version receives security fixes.

## Reporting a Vulnerability

**Please do not report security vulnerabilities through public GitHub Issues.**

Do not submit wet-lab validation data, proprietary sequences, unpublished constructs, patient data, confidential partner/customer data, private contact information, internal batch IDs, or exact confidential process parameters through public security reports or public GitHub Issues.

### Option 1 — GitHub Private Vulnerability Reporting (preferred)

Use GitHub's built-in private reporting:  
**Security → Report a vulnerability** on the [repository page](https://github.com/eijex/factorforge-cds/security/advisories/new).

This keeps the report confidential until a fix is released.

### Option 2 — Email

Send a report to **eijex.lab@gmail.com** with the subject line:  
`[SECURITY] FactorForge — <brief description>`

Include:
- FactorForge version affected (`pip show factorforge-cds`)
- Description of the vulnerability
- Steps to reproduce
- Potential impact
- Any suggested fix (optional)

### Response timeline

| Stage | Target time |
|-------|------------|
| Initial acknowledgement | 48 hours |
| Triage and severity assessment | 5 business days |
| Fix and patch release | Depends on severity |

We will notify you when the vulnerability is fixed and credit you in the release notes unless you prefer to remain anonymous.

## Scope

### In scope

- Input validation bypass in the API (`/api/optimize`)
- Sequence data exposure or logging
- Dependency vulnerabilities with exploitable attack vectors
- Remote code execution via crafted input sequences

### Out of scope

- **Dual-use research concerns** (e.g., "this tool could be used to design harmful proteins") — these are not security vulnerabilities in the software sense. FactorForge only optimizes the codon usage of a user-provided sequence; it does not design new biological functions. Dual-use concerns should be directed to the user's institutional biosafety committee.
- Performance issues or bugs without security impact
- Theoretical vulnerabilities without a proof of concept
- Issues in unsupported versions

## Security Design Notes

FactorForge is a stateless codon optimization tool:

- **No user data stored**: submitted sequences are not logged or persisted server-side
- **No authentication required**: the public API is read-only and stateless
- **Vercel serverless**: the web API runs on Vercel infrastructure; server-side security is managed by Vercel
- **Local use**: for sensitive or unpublished sequences, use the CLI (`pip install factorforge-cds`) or Docker image to run entirely offline
