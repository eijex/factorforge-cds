"""factorforge-validate CLI: convert wet-lab results into validation packages."""

import argparse
import hashlib
import sys
from pathlib import Path

from factorforge.validation import ValidationPackageGenerator
from factorforge.validation.package_generator import WetLabResult


def _sequence_hash(args: argparse.Namespace) -> str:
    if args.sequence_hash:
        if args.sequence_hash.startswith("sha256:"):
            return args.sequence_hash
        return f"sha256:{args.sequence_hash}"
    if args.sequence:
        return "sha256:" + hashlib.sha256(args.sequence.encode()).hexdigest()

    print("Error: --sequence or --sequence-hash is required.", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="factorforge-validate",
        description="Generate a structured validation package from wet-lab results.",
    )
    parser.add_argument(
        "--construct-id",
        required=True,
        help="Construct ID (e.g. CF-20260525-143022)",
    )
    parser.add_argument(
        "--version",
        required=True,
        help="FactorForge product version (e.g. 3.1.1)",
    )
    parser.add_argument("--host-profile", default="nbenthamiana")
    parser.add_argument("--profile", default=None)
    parser.add_argument("--sequence", default=None, help="CDS sequence (for hash only)")
    parser.add_argument(
        "--sequence-hash",
        default=None,
        help="Pre-computed sha256 hash (alternative to --sequence)",
    )
    parser.add_argument("--protein-name", required=True)
    parser.add_argument("--host-organism", default="N. benthamiana")
    parser.add_argument("--promoter", default=None)
    parser.add_argument("--subcellular-targeting", default=None)
    parser.add_argument("--expression-system", required=True)
    parser.add_argument("--harvest-timepoint", default=None)
    parser.add_argument(
        "--native-control",
        default=None,
        choices=["Yes", "No", "Not applicable"],
    )
    parser.add_argument(
        "--comparison",
        required=True,
        choices=["FactorForge better", "Equivalent", "Worse", "Not compared"],
    )
    parser.add_argument(
        "--expression-level",
        default=None,
        choices=["High", "Medium", "Low", "Not detected", "Not measured"],
    )
    parser.add_argument("--notes", default=None)
    parser.add_argument("--institution", default=None)
    parser.add_argument("--no-public", action="store_true", help="Opt out of public listing")
    parser.add_argument("--output-dir", default="validation_package", help="Output directory")

    args = parser.parse_args()
    result = WetLabResult(
        construct_id=args.construct_id,
        factorforge_version=args.version,
        host_profile=args.host_profile,
        profile=args.profile,
        sequence_hash=_sequence_hash(args),
        protein_name=args.protein_name,
        host_organism=args.host_organism,
        promoter=args.promoter,
        subcellular_targeting=args.subcellular_targeting,
        expression_system=args.expression_system,
        harvest_timepoint=args.harvest_timepoint,
        native_control=args.native_control,
        comparison=args.comparison,
        expression_level=args.expression_level,
        notes=args.notes,
        institution=args.institution,
        public_listing=not args.no_public,
    )

    output = ValidationPackageGenerator(Path(args.output_dir)).generate(result)
    print(f"Validation package generated: {output}")


if __name__ == "__main__":
    main()
