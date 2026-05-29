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
| `--gc-min` | `40` | Minimum GC% target |
| `--gc-max` | `55` | Maximum GC% target |
| `--format` | `fasta` | Output format: `fasta` or `genbank` |
| `--host` | `nbenthamiana` | Expression host: `nbenthamiana` or `by2` (Tobacco BY-2) |
| `--compare-profiles` | — | Comma-separated profiles to compare (e.g. `balanced,high_cai,gc_target`) |
| `--scan-mode` | `full` | Rule scan: `full` or `fast` |
| `--template` | — | MoClo construct template |

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
```
