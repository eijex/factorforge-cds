> **Note**: Internal workflow details are maintained in the private eijex-agentOps repository.

# FactorForge Agent Failure Modes

## Stop conditions

- Raw sequence detected in input or proposed output
- Forbidden claim phrase in generated text
- CI or test failure on an affected file
- Schema validation failure

## Degraded conditions

- Conflicting evidence on a biological parameter
- Ambiguous authorization state
- Schema drift between versions

## Recovery

Stop, report to a human, and await instruction. Never attempt autonomous
recovery after a claim or sequence violation.
