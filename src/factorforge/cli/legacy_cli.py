"""Archived v1 CLI compatibility shim."""

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Archived FactorForge v1 CLI shim. Use `factorforge optimize` instead."
    )
    parser.add_argument("args", nargs="*", help=argparse.SUPPRESS)
    parser.parse_args()
    raise SystemExit(
        "The v1 CLI is archived and is not shipped as a supported runtime. "
        "Use `factorforge optimize` for the supported profile/DP engines, or inspect "
        "`archive/v1-nbent-opticodon/` for historical provenance."
    )


if __name__ == "__main__":
    main()
