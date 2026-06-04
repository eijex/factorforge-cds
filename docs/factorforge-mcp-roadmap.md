# FactorForge MCP Roadmap

Model Context Protocol servers are planned, not implemented. They will connect
AI agents such as Claude Code and IDE assistants to FactorForge design tools and
biological evidence sources through the MCP open standard.

Every server listed here is a roadmap item. No MCP server should be documented
as available until it has implementation, tests, and user-facing setup docs.

## Planned Servers

### factorforge-sequence-mcp

Status: planned, not implemented.

Planned tools:

- `optimize_sequence`
- `reverse_translate`
- `score_sequence`
- `compare_sequence_candidates`
- `export_fasta`
- `export_genbank`

This server will expose sequence-design operations for controlled agent use.
Tool outputs should include profile, host, score components, rule findings, and
warnings about unavailable optional calculations.

### factorforge-host-profile-mcp

Status: planned, not implemented.

Planned tools:

- `list_host_profiles`
- `get_host_profile`
- `validate_host_profile`
- `compare_host_profiles`
- `create_profile_draft`

This server will expose host-profile registry operations. It should enforce the
profile schema and keep validation status separate from profile availability.

### factorforge-rule-audit-mcp

Status: planned, not implemented.

Planned tools:

- `scan_homopolymers`
- `scan_gc_windows`
- `scan_restriction_sites`
- `scan_polyA_motifs`
- `scan_rare_codon_runs`
- `generate_rule_audit_report`

This server will expose sequence rule-audit workflows. Reports should include
thresholds, scan windows, motif lists, and severity semantics.

### factorforge-evidence-mcp

Status: planned, not implemented.

Planned tools:

- `pubmed_search`
- `ncbi_gene_lookup`
- `ncbi_nucleotide_lookup`
- `build_host_profile_evidence_table`
- `build_codon_usage_source_audit`

This server will connect profile design to evidence review. It must distinguish
source audit, literature notes, and validation claims.

### factorforge-dev-inbox-mcp

Status: planned, not implemented.

Planned tools:

- `create_github_issue`
- `summarize_ci_failure`
- `prepare_claude_handoff`
- `generate_pr_review_checklist`

This server will support developer workflow packaging. It should prepare context
for review, not automatically fix, merge, deploy, or change biological scoring
semantics.
