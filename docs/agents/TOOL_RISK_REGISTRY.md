> **Note**: Internal workflow details are maintained in the private eijex-agentOps repository.

# FactorForge Agent Tool Risk Registry

| Tool category | Risk level | Notes |
|---|---|---|
| Public metadata query | low | Read-only external query |
| Public result read | low | Reads supported public output |
| Design computation | medium | Requires human confirmation of inputs |
| Artifact write | high | Requires explicit human confirmation |
| Public publishing action | critical | Requires explicit human authorization |

Risk levels are low, medium, high, and critical. Five consecutive high or
critical actions require renewed human confirmation.
