#!/usr/bin/env python3
"""Deterministic sfGFP worked example for FactorForge.

Usage:
    python run_example.py           # verify frozen output matches current run
    python run_example.py --freeze  # overwrite frozen output (first run only)

Profile: assembly_friendly | Seed: 320 | Host: N. benthamiana
Scoring contract: v1.1 (multi_constraint_pass = biological_pass AND assembly_pass AND gc_in_target_range)
Reference: Pédelacq et al. 2006, Nat Biotechnol 24:79-88, PMID 16369541 | PDB:2B3P
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import factorforge as _ff
from factorforge.engines.profile.optimizer import RuleBasedOptimizer
from factorforge.engines.profile.rules.domesticator import Domesticator
from factorforge.utils.sequence_validator import validate_cds_output

PROFILE = "assembly_friendly"
SEED = 320
GC_MIN = 55.0  # N. benthamiana benchmark config (scoring_contract v1.1)
GC_MAX = 65.0

_DOM = Domesticator()


def load_sfgfp_aa() -> str:
    faa = HERE / "input_sequence.faa"
    lines = faa.read_text(encoding="utf-8").strip().split("\n")
    return "".join(l for l in lines if not l.startswith(">")).replace(" ", "")


def compute_validation_metrics(aa_seq: str, dna_seq: str, gc_percent: float) -> dict:
    """Derive scoring_contract v1.1 primitives from production factorforge functions.

    NOTE: OptimizationResult.metadata does NOT contain biological_pass or assembly_pass.
    These must be derived explicitly using production validators.
    """
    val = validate_cds_output(aa_seq, dna_seq)
    biological_pass = bool(val["passed"])

    type_iis = _DOM.scan_restriction_sites(dna_seq, "golden_gate")
    assembly_pass = len(type_iis) == 0
    forbidden_count = len(type_iis)

    gc_in_target_range = GC_MIN <= gc_percent <= GC_MAX

    # multi_constraint_pass is DERIVED — never assert directly
    multi_constraint_pass = biological_pass and assembly_pass and gc_in_target_range

    return {
        "biological_pass": biological_pass,
        "assembly_pass": assembly_pass,
        "forbidden_type_iis_site_count": forbidden_count,
        "gc_in_target_range": gc_in_target_range,
        "multi_constraint_pass": multi_constraint_pass,
        "aa_identity": round(float(val.get("aa_identity", 1.0)), 4),
    }


def build_design_package(result, aa_seq: str) -> dict:
    """Build design_package.json conforming to design_package.schema.json v1.0."""
    seq_hash = hashlib.sha256(result.sequence.encode()).hexdigest()
    gc = result.metrics.get("gc_percent", 0.0)
    cai = result.metrics.get("cai", 0.0)
    v = compute_validation_metrics(aa_seq, result.sequence, gc)
    return {
        "design_package_version": "1.0",
        "design_id": f"sfgfp_{PROFILE}_seed{SEED}_v1",
        "sequence_type": "cds",
        "input_type": "protein",
        "host_profile": {
            "id": "nbenthamiana",
            "display_name": "Tobacco (N. benthamiana)",
            "scientific_name": "Nicotiana benthamiana",
            "ncbi_taxonomy_id": "4100",
            "status": "stable",
        },
        "optimization": {"engine": "profile", "profile": PROFILE},
        "metrics": {
            "cai": round(cai, 4),
            "gc_percent": round(gc, 2),
            "mfe_kcal_mol": None,
            "mfe_status": "not_computed",
            "mfe_used": False,
        },
        "validation": {
            "biological_pass": v["biological_pass"],
            "assembly_pass": v["assembly_pass"],
            "aa_identity": v["aa_identity"],
            "forbidden_type_iis_site_count": v["forbidden_type_iis_site_count"],
        },
        "evidence": {
            "sequence_hash": seq_hash,
            "registry_version": getattr(_ff, "__version__", "3.2.0"),
            "registry_hash": "N/A",
        },
        "claim_boundary": {
            "in_silico_only": True,
            "no_yield_claim": True,
            "no_wet_lab_claim": True,
            "no_clinical_claim": True,
        },
    }


def build_validation_summary(design_pkg: dict) -> dict:
    """Build validation_summary.json for ValidationHub computational intake.

    Reads primitives from design_pkg to avoid recomputation.
    multi_constraint_pass is re-derived here as canonical check.
    """
    bio_pass = design_pkg["validation"]["biological_pass"]
    asm_pass = design_pkg["validation"]["assembly_pass"]
    gc = design_pkg["metrics"]["gc_percent"]
    gc_in_target_range = GC_MIN <= gc <= GC_MAX
    multi_constraint_pass = bio_pass and asm_pass and gc_in_target_range  # derived
    return {
        "schema_version": "0.1.0",
        "design_package_id": design_pkg["design_id"],
        "sequence_publication_source": (
            "Pédelacq et al. 2006, Nat Biotechnol 24:79-88, PMID 16369541"
        ),
        "public_example_id": "sfgfp_assembly_friendly_seed320_nbenthamiana_v1",
        "host_organism": "Nicotiana benthamiana",
        "evidence_type": "computational_validation",
        "computational": {
            "scoring_contract_version": "v1.1",
            "gc_in_target_range": gc_in_target_range,
            "multi_constraint_pass": multi_constraint_pass,
            "multi_constraint_pass_definition": (
                "biological_pass AND assembly_pass AND gc_in_target_range"
            ),
        },
        "experimental": {
            "assay_type": "not_provided",
            "outcome_status": "not_tested",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--freeze",
        action="store_true",
        help="Overwrite frozen output files (first run only)",
    )
    args = parser.parse_args()

    aa_seq = load_sfgfp_aa()
    optimizer = RuleBasedOptimizer()
    result = optimizer.optimize(aa_seq, profile=PROFILE, seed=SEED)

    design_pkg = build_design_package(result, aa_seq)
    val_summary = build_validation_summary(design_pkg)

    out_dir = HERE / "output"
    out_dir.mkdir(exist_ok=True)
    dp_path = out_dir / "design_package.json"
    vs_path = out_dir / "validation_summary.json"

    if args.freeze:
        dp_path.write_text(
            json.dumps(design_pkg, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        vs_path.write_text(
            json.dumps(val_summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"Frozen: {dp_path}")
        print(f"Frozen: {vs_path}")
        print(f"\nKey results:")
        print(f"  CAI:                 {design_pkg['metrics']['cai']}")
        print(f"  GC%:                 {design_pkg['metrics']['gc_percent']}")
        print(f"  gc_in_target_range:  {val_summary['computational']['gc_in_target_range']}")
        print(f"  multi_constraint_pass: {val_summary['computational']['multi_constraint_pass']}")
    else:
        errors = []
        for path, obj in [(dp_path, design_pkg), (vs_path, val_summary)]:
            if not path.exists():
                errors.append(f"Missing frozen output: {path}")
                continue
            frozen = json.loads(path.read_text(encoding="utf-8"))
            if frozen != obj:
                errors.append(f"Reproducibility drift in {path.name}")
        if errors:
            for e in errors:
                print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        else:
            print("OK: frozen outputs match current run (reproducible)")


if __name__ == "__main__":
    main()
