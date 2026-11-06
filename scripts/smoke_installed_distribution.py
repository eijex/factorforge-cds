"""Check the public Bioconda distribution, never an editable checkout.

Run with the installed environment's ``python -I``. Only synthetic inputs are used.
"""

import argparse
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from urllib.parse import urlparse


def verify_origin(module_path, prefix, version):
    """Fail if imports escape the prefix or conda provenance does not match."""
    prefix = Path(prefix).resolve()
    module_path = Path(module_path).resolve()
    if not module_path.is_relative_to(prefix):
        raise ValueError("Imported module is outside the installed environment")
    records = list((prefix / "conda-meta").glob("factorforge-cds-*.json"))
    if len(records) != 1:
        raise ValueError("Expected exactly one installed conda package record")
    record = json.loads(records[0].read_text(encoding="utf-8"))
    if record.get("name") != "factorforge-cds" or record.get("version") != version:
        raise ValueError("Installed conda record does not match imported version")
    channel = record.get("channel", "")
    url = urlparse(channel)
    if not (
        channel == "bioconda"
        or (
            url.scheme == "https"
            and url.hostname in {"conda.anaconda.org", "conda.anaconda.com"}
            and url.path.rstrip("/") in {"/bioconda", "/bioconda/noarch"}
        )
    ):
        raise ValueError("Installed conda record is not from the Bioconda channel")
    # Confirm this module actually belongs to the installed artifact's file list.
    relative = module_path.relative_to(prefix).as_posix()
    if relative not in record.get("files", []):
        raise ValueError("Imported module is not listed in the installed conda artifact")
    return {key: record.get(key) for key in ("name", "version", "build", "channel", "sha256")}


def exercise(report):
    from Bio import SeqIO
    from Bio.Seq import Seq

    import factorforge

    version = importlib.metadata.version("factorforge-cds")
    report.update(version=version, python=sys.version, module_path=factorforge.__file__)
    report["package"] = verify_origin(factorforge.__file__, sys.prefix, version)
    executable = Path(sys.prefix) / "bin" / "factorforge"
    if os.name == "nt":
        executable = Path(sys.prefix) / "Scripts" / "factorforge.exe"
    if not executable.is_file():
        raise ValueError(f"Installed CLI launcher is missing: {executable}")
    report["cli_path"] = str(executable)
    with tempfile.TemporaryDirectory(prefix="factorforge-distribution-") as temporary:
        cwd = Path(temporary)

        def run(*args):
            result = subprocess.run(
                [str(executable), *args],
                cwd=cwd,
                text=True,
                encoding="utf-8",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=120,
                check=False,
                env={
                    key: value
                    for key, value in os.environ.items()
                    if key.upper() not in {"PYTHONPATH", "PYTHONHOME"}
                },
            )
            report["calls"].append(
                {"args": args, "exit_code": result.returncode, "output": result.stdout}
            )
            return result

        for args in [("--version",), ("--help",), ("optimize", "--help")]:
            result = run(*args)
            if result.returncode or (args == ("--version",) and version not in result.stdout):
                raise ValueError("Installed CLI version/help failed; inspect calls")
        protein = "MASTWQHLEK" * 6
        (cwd / "input.fasta").write_text(">synthetic_smoke\n" + protein + "\n", encoding="utf-8")
        sequences = []
        for name, options in [("default.fasta", []), ("explicit.fasta", ["--engine", "dp"])]:
            result = run("optimize", "input.fasta", *options, "-o", name)
            if result.returncode:
                raise ValueError("Installed DP CLI failed; inspect calls")
            records = list(SeqIO.parse(cwd / name, "fasta"))
            if len(records) != 1 or "engine=dp" not in records[0].description:
                raise ValueError("Expected exactly one DP output record")
            dna = str(records[0].seq).upper()
            if not set(dna) <= set("ACGT") or len(dna) % 3:
                raise ValueError("Invalid CDS alphabet or reading frame")
            if str(Seq(dna).translate()).removesuffix("*") != protein:
                raise ValueError("Amino-acid identity failed")
            sequences.append(dna)
        if sequences[0] != sequences[1]:
            raise ValueError("Default and explicit DP output differ")
        (cwd / "invalid.fasta").write_text(">invalid\nMAST!XYZ\n", encoding="utf-8")
        if run("optimize", "invalid.fasta", "-o", "invalid-output.fasta").returncode == 0:
            raise ValueError("Invalid protein input was accepted")
        report.update(
            translation_identity=True,
            default_matches_explicit_dp=True,
            invalid_input_rejected=True,
            protein_length=len(protein),
            output_nt=len(sequences[0]),
            output_sha256=hashlib.sha256(sequences[0].encode()).hexdigest(),
        )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "FAIL",
        "calls": [],
        "scope": "Installed DP smoke only; not 3DP/4DP, live web or biological validation",
    }
    try:
        exercise(report)
        report["status"] = "PASS"
    except (ValueError, OSError, ImportError, subprocess.SubprocessError) as error:
        report["error"] = f"{type(error).__name__}: {error}"
    (args.output / "smoke-result.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "calls"}, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
