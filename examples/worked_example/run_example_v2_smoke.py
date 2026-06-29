#!/usr/bin/env python3
"""Current-default sfGFP smoke check for FactorForge.

Job 168 / v3.3.0 (_analysis/025) introduced a new production-default codon
reference (NbeV1.1 LAB-strain, GC 40-47%) for N. benthamiana, replacing the
legacy reference that run_example.py is pinned to. That promotion was
provisionally reverted on 2026-06-29 pending an MFE re-sensitivity + 2x2
factorial recheck (see data/reference/active_codon_reference.json) — this
script does not hardcode which one is active; it tracks whatever the engine
currently resolves to.

This script does NOT compare against a frozen artifact — it only verifies
that the current production default runs successfully end-to-end and
reports the provenance facts (resolved asset ID, contract version, GC band)
a reader needs to interpret the output correctly. Promoting this to a
frozen, version-pinned worked example is a separate, deliberate follow-up
decision — not a byproduct of this job.

Usage:
    python run_example_v2_smoke.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from factorforge.engines.profile.optimizer import RuleBasedOptimizer  # noqa: E402
from factorforge.utils.sequence_validator import validate_cds_output  # noqa: E402

PROFILE = "assembly_friendly"
SEED = 320
ACTIVE_REFERENCE_PATH = ROOT / "data" / "reference" / "active_codon_reference.json"


def load_sfgfp_aa() -> str:
    faa = HERE / "input_sequence.faa"
    lines = faa.read_text(encoding="utf-8").strip().split("\n")
    return "".join(line for line in lines if not line.startswith(">")).replace(" ", "")


def main() -> None:
    aa_seq = load_sfgfp_aa()

    # No codon_table_path override ⇒ tracks whatever the engine's current
    # production default actually is (see engines/profile/utils.py
    # _HOST_CODON_TABLE_OVERRIDES).
    optimizer = RuleBasedOptimizer()
    result = optimizer.optimize(aa_seq, profile=PROFILE, seed=SEED)

    val = validate_cds_output(aa_seq, result.sequence)
    if not val["passed"]:
        print(f"ERROR: biological validation failed: {val}", file=sys.stderr)
        sys.exit(1)

    active_ref = json.loads(ACTIVE_REFERENCE_PATH.read_text(encoding="utf-8"))
    contract_version = active_ref["codon_reference_contract_version"]

    print(f"OK: current-default ({contract_version}) run completed successfully.")
    print("\nProvenance:")
    print(f"  active_codon_table_id:        {active_ref['active_codon_table_id']}")
    print(f"  active_asset_type:            {active_ref['active_asset_type']}")
    print(f"  codon_reference_contract_version: {active_ref['codon_reference_contract_version']}")
    print(f"  gc_range:                     {active_ref['gc_range']}")
    print("\nResult:")
    print(f"  CAI:                 {result.metrics.get('cai')}")
    print(f"  GC%:                 {result.metrics.get('gc_percent')}")
    print(f"  aa_identity:          {val.get('aa_identity', 1.0)}")
    print(f"  sequence length (nt): {len(result.sequence)}")


if __name__ == "__main__":
    main()
