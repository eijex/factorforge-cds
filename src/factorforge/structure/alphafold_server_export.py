"""Offline AlphaFold Server (alphafoldserver.com) job JSON exporter.

Generates validated multi-chain job JSON matching Google DeepMind's AlphaFold
Server web schema without making any network calls or requiring local GPU.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Offline AlphaFold Server export · No network calls made
"""

from __future__ import annotations

import json
import re
from collections.abc import Sequence
from pathlib import Path, PureWindowsPath
from tempfile import NamedTemporaryFile
from typing import Any

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
MIN_PROTEIN_LENGTH = 4
MAX_JOB_TOKENS = 5000
MAX_MODEL_SEED = 2**32 - 1


class AlphaFoldExportError(ValueError):
    """Raised when sequence validation or AlphaFold Server schema constraints fail."""


def sanitize_job_name(name: str | None, max_length: int = 64) -> str:
    """Sanitize a job name to safe alphanumeric and underscore/hyphen characters.

    Strips directory components and control characters from the supplied name.
    The caller is responsible for choosing a non-confidential name.
    """
    if not name:
        return "factorforge_job"
    # Take basename if a path was inadvertently supplied
    base = PureWindowsPath(name).name
    # Keep only safe characters: letters, digits, underscores, hyphens
    sanitized = re.sub(r"[^A-Za-z0-9_-]", "_", base)
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    if not sanitized:
        return "factorforge_job"
    return sanitized[:max_length]


def validate_protein_sequence(sequence: str, chain_label: str = "Chain") -> str:
    """Validate and normalize a protein sequence for AlphaFold Server intake.

    Args:
        sequence: Raw amino acid sequence string.
        chain_label: Identifier for error reporting.

    Returns:
        Normalized uppercase amino acid sequence string.

    Raises:
        AlphaFoldExportError: If sequence is empty, out of bounds, or contains
            non-standard amino acids.
    """
    if not isinstance(sequence, str):
        raise AlphaFoldExportError(f"{chain_label} sequence must be text.")
    residues = "".join(sequence.split())
    if not residues.isascii():
        raise AlphaFoldExportError(f"{chain_label} contains an invalid or non-standard amino acid.")
    clean_seq = residues.upper()
    if not clean_seq:
        raise AlphaFoldExportError(f"{chain_label} sequence is empty.")

    if len(clean_seq) < MIN_PROTEIN_LENGTH:
        raise AlphaFoldExportError(
            f"{chain_label} sequence length ({len(clean_seq)} aa) is below minimum {MIN_PROTEIN_LENGTH} aa."
        )

    if len(clean_seq) > MAX_JOB_TOKENS:
        raise AlphaFoldExportError(
            f"{chain_label} sequence length ({len(clean_seq)} aa) exceeds maximum {MAX_JOB_TOKENS} aa."
        )

    if set(clean_seq) - STANDARD_AMINO_ACIDS:
        raise AlphaFoldExportError(
            f"{chain_label} contains an invalid or non-standard amino acid. "
            "AlphaFold Server requires standard 20 amino acids (ACDEFGHIKLMNPQRSTVWY)."
        )

    return clean_seq


def parse_protein_fasta(content_or_path: str | Path) -> list[str]:
    """Parse one or more protein sequences from a FASTA file or string.

    Omits FASTA headers and metadata from the returned sequences to prevent
    information leakage.

    Returns:
        List of normalized amino acid sequences.
    """
    # A Path explicitly requests file IO; a string is FASTA content.
    if isinstance(content_or_path, Path):
        text = content_or_path.read_text(encoding="utf-8-sig")
    else:
        text = content_or_path.removeprefix("\ufeff")

    sequences: list[str] = []
    current_seq_parts: list[str] = []
    seen_header = False

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seen_header:
                sequences.append("".join(current_seq_parts))
                current_seq_parts = []
            seen_header = True
        else:
            if not seen_header:
                raise AlphaFoldExportError("Expected a FASTA header before sequence data.")
            current_seq_parts.append(line)

    if seen_header:
        sequences.append("".join(current_seq_parts))

    if not sequences:
        raise AlphaFoldExportError("No FASTA sequence records found in input.")

    return [validate_protein_sequence(s, f"FASTA record {i + 1}") for i, s in enumerate(sequences)]


def build_alphafold_server_job_payload(
    chains: Sequence[tuple[str, int] | dict[str, Any] | str],
    job_name: str | None = None,
    model_seeds: Sequence[int] | int | None = 1,
) -> list[dict[str, Any]]:
    """Build a validated AlphaFold Server job payload matching the server schema.

    Schema:
    [
      {
        "name": "<job_name>",
        "dialect": "alphafoldserver",
        "version": 1,
        "modelSeeds": ["1"],
        "sequences": [
          {
            "proteinChain": {
              "sequence": "<AMINO_ACIDS>",
              "count": 1
            }
          }
        ]
      }
    ]

    Args:
        chains: Sequence of chains, where each item can be:
            - tuple of (sequence, count)
            - dict with keys 'sequence' and optional 'count' (default 1)
            - string containing the sequence (count defaults to 1)
        job_name: Optional job name string.
        model_seeds: uint32 seed or sequence of seeds (default: 1).
            Pass None or an empty sequence for server-assigned seeds.

    Returns:
        List containing the single formatted job dictionary.
    """
    if not chains:
        raise AlphaFoldExportError("At least one protein chain is required.")

    sanitized_name = sanitize_job_name(job_name)

    # Normalize seeds
    if isinstance(model_seeds, int):
        seeds = [model_seeds]
    elif model_seeds is None:
        seeds = []
    elif isinstance(model_seeds, Sequence) and not isinstance(model_seeds, (str, bytes)):
        seeds = list(model_seeds)
    else:
        raise AlphaFoldExportError("Expected an integer seed or a sequence of integer seeds.")
    if any(type(seed) is not int or not 0 <= seed <= MAX_MODEL_SEED for seed in seeds):
        raise AlphaFoldExportError(f"Model seeds must be integers between 0 and {MAX_MODEL_SEED}.")

    sequences_list: list[dict[str, Any]] = []
    total_tokens = 0

    for idx, item in enumerate(chains, start=1):
        chain_seq: str
        chain_count: int = 1

        if isinstance(item, tuple) and len(item) == 2:
            chain_seq, chain_count = item
        elif isinstance(item, dict):
            if "sequence" not in item:
                raise AlphaFoldExportError(f"Chain {idx} dict missing required 'sequence' key.")
            chain_seq = item["sequence"]
            chain_count = item.get("count", 1)
        elif isinstance(item, str):
            chain_seq = item
            chain_count = 1
        else:
            raise AlphaFoldExportError(f"Unsupported chain descriptor type at index {idx}: {type(item)}")

        if type(chain_count) is not int or chain_count < 1:
            raise AlphaFoldExportError(f"Chain {idx} count must be >= 1 and an integer.")

        valid_seq = validate_protein_sequence(chain_seq, f"Chain {idx}")
        total_tokens += len(valid_seq) * chain_count
        if total_tokens > MAX_JOB_TOKENS:
            raise AlphaFoldExportError(f"Total protein tokens including copies exceed {MAX_JOB_TOKENS}.")

        sequences_list.append({
            "proteinChain": {
                "sequence": valid_seq,
                "count": chain_count,
            }
        })

    job_dict: dict[str, Any] = {
        "name": sanitized_name,
        "dialect": "alphafoldserver",
        "version": 1,
        "modelSeeds": [str(seed) for seed in seeds],
        "sequences": sequences_list,
    }

    return [job_dict]


def export_alphafold_server_job(
    chains: Sequence[tuple[str, int] | dict[str, Any] | str],
    output_path: Path | str | None = None,
    job_name: str | None = None,
    model_seeds: Sequence[int] | int | None = 1,
    overwrite: bool = False,
) -> tuple[list[dict[str, Any]], str]:
    """Build and optionally export an AlphaFold Server job JSON file.

    Args:
        chains: Protein chain definitions.
        output_path: Optional file path to write.
        job_name: Name for the job manifest.
        model_seeds: Random seed(s).
        overwrite: Whether to overwrite existing destination file.

    Returns:
        Tuple of (parsed payload list, formatted JSON string).

    Raises:
        FileExistsError: If output_path exists and overwrite=False.
        AlphaFoldExportError: On validation failure.
    """
    payload = build_alphafold_server_job_payload(
        chains=chains,
        job_name=job_name,
        model_seeds=model_seeds,
    )
    json_str = json.dumps(payload, indent=2) + "\n"

    if output_path is not None:
        target = Path(output_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        if overwrite:
            with NamedTemporaryFile(dir=target.parent, mode="w", encoding="utf-8", delete=False) as stream:
                temporary = Path(stream.name)
                try:
                    stream.write(json_str)
                    stream.close()
                    temporary.replace(target)
                finally:
                    stream.close()
                    temporary.unlink(missing_ok=True)
        else:
            try:
                with target.open("x", encoding="utf-8") as stream:
                    stream.write(json_str)
            except FileExistsError as exc:
                raise FileExistsError("Output already exists. Use --overwrite to replace it.") from exc

    return payload, json_str
