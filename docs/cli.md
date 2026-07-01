# CLI Reference

## Commands

### `factorforge optimize`

Optimize a protein sequence into a host-compatible CDS (*N. benthamiana* by default; `--host by2` for Tobacco BY-2).

```bash
factorforge optimize input.fasta -o output.fasta
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--engine`, `-e` | `dp` | Engine: `dp` (feasibility) or `profile` (rule/profile) |
| `--profile`, `-p` | `balanced` | Optimization profile |
| `--objective` | `feasibility_best` | DP objective |
| `--gc-min` | `55` | Minimum GC% target |
| `--gc-max` | `65` | Maximum GC% target |
| `--format` | `fasta` | Output format: `fasta` or `genbank` |
| `--host` | `nbenthamiana` | Expression host: `nbenthamiana` or `by2` (Tobacco BY-2) |
| `--reference-id` | — | Expert/research codon-reference ID; checksum-validated and not recommended for production defaults |
| `--compare-profiles` | — | Comma-separated profiles to compare (e.g. `balanced,high_cai,gc_target`) |
| `--scan-mode` | `full` | Rule scan: `full` or `fast` |
| `--template` | — | MoClo construct template |

`--reference-id` is an expert/research path for replaying or comparing packaged
codon-reference tables from `data/reference/reference_policy_manifest.json`.
Each selected table is checksum-validated before use. Non-production tiers
continue after printing a warning that quotes the manifest claim boundary; they
are not production recommendations. Public REST API reference overrides remain
unsupported.

## API Endpoints

The web API (`https://factorforge.eijex.com`) exposes the following endpoints:

| Endpoint | Description |
|----------|-------------|
| `POST /api/optimize` | Single sequence optimization |
| `POST /api/optimize/compare` | Compare multiple profiles side-by-side |
| `POST /api/optimize/batch` | Optimize up to 20 sequences in one request |

For AI agent access, use [Eijex MCP](https://mcp.eijex.com) which wraps these endpoints as MCP tools.

### `factorforge list-engines`

List all available optimization engines.

```bash
factorforge list-engines
```

## Examples

```bash
# Default (DP feasibility engine)
factorforge optimize input.fasta -o output.fasta

# Profile engine
factorforge optimize input.fasta -e profile -p balanced -o output.fasta

# GenBank output with MoClo template
factorforge optimize input.fasta -e profile -p balanced \
  --template standard_expression -o output.gb --format genbank

# Custom GC target range
factorforge optimize input.fasta --gc-min 45 --gc-max 60 -o output.fasta

# Fast scan (skip rare codon run detection)
factorforge optimize input.fasta --scan-mode fast -o output.fasta

# Expert/research reference replay or comparison
factorforge optimize input.fasta --reference-id nbenthamiana_qld183_v103
```
