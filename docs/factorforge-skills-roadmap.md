# FactorForge Skills Roadmap

Agent skills are planned, not implemented. They should package repeatable
review workflows for AI coding agents while keeping sequence-design decisions
explicit and auditable.

No skill listed here should be treated as available until it has implementation,
tests or dry-run examples, and user-facing usage documentation.

## Planned Skills

### host-profile-design-skill

Status: planned, not implemented.

Purpose: scaffold a new host profile YAML with a schema validation checklist.

Expected scope:

- Request species, strain or cell line, tissue or cell type, and expression
  platform separately.
- Require codon usage source and table status.
- Require validation status and not-recommended-for notes.
- Prevent profile drafts from implying wet-lab validation.

### sequence-calibration-skill

Status: planned, not implemented.

Purpose: provide a checklist for safely reviewing CAI, GC, homopolymer, and
repair-threshold changes.

Expected scope:

- Identify changed constants and profile values.
- Require fixture-based score comparison.
- Require documentation updates when scoring semantics change.
- Flag user-visible behavior changes before release.

### assembly-friendly-skill

Status: planned, not implemented.

Purpose: review Golden Gate, Type IIS, and synthesis-friendly design scope.

Expected scope:

- Confirm which restriction enzymes are scanned.
- Confirm whether enzyme lists are configurable.
- Check CAI guardrails after synonymous repair.
- Separate implemented constraints from planned constraints.

### rule-audit-skill

Status: planned, not implemented.

Purpose: check threshold inconsistencies, magic numbers, severity semantics, and
rule-reporting clarity.

Expected scope:

- Compare rule constants with documentation.
- Verify rule findings include thresholds and context.
- Detect duplicated thresholds across modules.
- Review false-positive notes for heuristic motif scanners.

### factorforge-dev-agent-skill

Status: planned, not implemented.

Purpose: convert GitHub Issues and developer-inbox items into Claude Code
implementation prompts.

Expected scope:

- Summarize observed failure or request.
- Identify affected files and tests.
- Preserve explicit safety boundaries.
- Avoid auto-fix, auto-merge, auto-deploy, or unreviewed scoring changes.

### documentation-skill

Status: planned, not implemented.

Purpose: ensure docs stay synchronized after code changes involving scoring,
host profiles, rule thresholds, and workflow boundaries.

Expected scope:

- Check docs against live constants.
- Confirm roadmap items are marked as planned when not implemented.
- Update examples when profile schemas change.
- Flag missing docs for user-visible behavior changes.
