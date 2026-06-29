"""
Utility helpers for FactorForge profile engine engines.
"""

from __future__ import annotations

import importlib.resources
import json
import os
from pathlib import Path
from typing import Any, cast


def get_data_path() -> Path:
    """Get data directory path from environment or package resources.

    Checks FACTORFORGE_DATA_DIR environment variable first, then falls back to
    the bundled package data directory (works with pip-installed packages).

    Returns:
        Path to data directory.

    Examples:
        >>> data_dir = get_data_path()
        >>> isinstance(data_dir, Path)
        True
    """
    data_dir = os.getenv("FACTORFORGE_DATA_DIR")
    if data_dir:
        return Path(data_dir)

    # Use importlib.resources to locate bundled package data
    # This works correctly whether the package is installed via pip or run from source
    try:
        ref = importlib.resources.files("factorforge") / "data"
        return Path(str(ref))
    except (TypeError, ModuleNotFoundError):
        # Fallback for development: go up from profile/utils.py to src/factorforge/data
        return Path(__file__).resolve().parents[2] / "data"


def calculate_gc(sequence: str) -> float:
    """Calculate GC content percentage.

    Args:
        sequence: DNA sequence string.

    Returns:
        GC content as percentage (0-100).

    Examples:
        >>> calculate_gc("ATGC")
        50.0
    """
    if not sequence:
        return 0.0

    seq = sequence.upper()
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(sequence)) * 100


def count_dinucleotides(sequence: str, dinucleotide: str = "CG") -> int:
    """Count occurrences of a dinucleotide in a DNA sequence.

    Args:
        sequence: DNA sequence string (case-insensitive).
        dinucleotide: Two-character dinucleotide to count (e.g., "CG", "TA").

    Returns:
        Count of dinucleotide occurrences.

    Examples:
        >>> count_dinucleotides("ACGACG", "CG")
        2
    """
    seq = sequence.upper()
    dn = dinucleotide.upper()
    return sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == dn)


def calculate_dinucleotide_ratio(sequence: str, dinucleotide: str = "CG") -> float:
    """Calculate observed/expected ratio of a dinucleotide.

    Compares actual dinucleotide frequency to what would be expected
    from mononucleotide composition. Ratio < 1.0 means suppressed,
    > 1.0 means enriched.

    Args:
        sequence: DNA sequence string (case-insensitive).
        dinucleotide: Two-character dinucleotide (e.g., "CG", "TA").

    Returns:
        Observed/expected ratio (0.0 if sequence too short or denominator zero).

    Examples:
        >>> calculate_dinucleotide_ratio("ACGTACGT", "CG")  # doctest: +SKIP
        1.0
    """
    seq = sequence.upper()
    if len(seq) < 2:
        return 0.0

    dn = dinucleotide.upper()
    observed = sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == dn)

    n1 = seq.count(dn[0])
    n2 = seq.count(dn[1])
    n = len(seq)

    expected = (n1 * n2) / n if n > 0 else 0.0
    if expected == 0.0:
        return 0.0

    return observed / expected


# Job 168 / v3.3.0 (_analysis/025) introduced a host -> production-default
# codon table file override mechanism and pointed nbenthamiana at the NbeV1.1
# LAB-strain derived table (released as part of v3.2.7). Provisionally
# reverted to empty (falls back to the legacy {host}_codons.json convention)
# pending an MFE re-sensitivity + 2x2 factorial recheck. The NbeV1.1 table
# remains on disk and selectable; see data/reference/active_codon_reference.json.
_HOST_CODON_TABLE_OVERRIDES: dict[str, str] = {}


def resolve_host_codon_table_path(host: str, codon_tables_dir: Path) -> Path:
    """Resolve the production-default codon table file path for a host."""
    override = _HOST_CODON_TABLE_OVERRIDES.get(host)
    filename = override or f"{host}_codons.json"
    return codon_tables_dir / filename


def load_codon_table(organism: str, codon_tables_dir: Path) -> dict[str, Any]:
    """Load codon usage table for organism.

    Args:
        organism: Organism name (e.g., "human", "ecoli").
        codon_tables_dir: Directory containing codon table files.

    Returns:
        Codon table payload parsed from JSON.

    Raises:
        FileNotFoundError: If codon table file not found.
    """
    if organism.endswith(".json"):
        codon_table_path = codon_tables_dir / organism
    else:
        codon_table_path = resolve_host_codon_table_path(organism, codon_tables_dir)

    with open(codon_table_path, "r", encoding="utf-8") as handle:
        return cast(dict[str, Any], json.load(handle))


def build_aa_to_codons_map(codon_table: dict[str, Any]) -> dict[str, list[str]]:
    """Build amino-acid-to-codons map from a codon table payload."""
    amino_acids = codon_table.get("amino_acids", {})
    if not isinstance(amino_acids, dict):
        return {}
    return {aa: list(info.get("codons", [])) for aa, info in amino_acids.items()}


def load_golden_set(data_dir: Path | None = None) -> dict[str, Any]:
    """Load golden set codon table for CAI reference weights.

    Falls back to the standard codon table if golden set file is not found.

    Args:
        data_dir: Data directory path. Defaults to get_data_path().

    Returns:
        Golden set codon table dict.
    """
    if data_dir is None:
        data_dir = get_data_path()
    golden_path = data_dir / "nbenthamiana_golden_set.json"
    if golden_path.exists():
        with open(golden_path, "r", encoding="utf-8") as f:
            return cast(dict[str, Any], json.load(f))
    # Fallback to standard table
    standard_path = data_dir / "nbenthamiana_codons.json"
    with open(standard_path, "r", encoding="utf-8") as f:
        return cast(dict[str, Any], json.load(f))


def translate_codon(codon: str, codon_table: dict[str, str]) -> str:
    """Translate DNA codon to amino acid.

    Args:
        codon: 3-letter DNA codon.
        codon_table: Codon to amino acid mapping.

    Returns:
        Single letter amino acid code.

    Raises:
        KeyError: If codon not in table.
    """
    return codon_table[codon.upper()]


def parse_fasta_records(content: str) -> list[tuple[str, str]]:
    """Parse FASTA content into (record_id, sequence) tuples.

    Args:
        content: FASTA text content.

    Returns:
        List of (record_id, sequence) tuples.

    Raises:
        ValueError: If content is not valid FASTA.
    """
    records: list[tuple[str, str]] = []
    seq_id: str | None = None
    seq_lines: list[str] = []

    for raw_line in content.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if seq_id is not None:
                sequence = "".join(seq_lines).upper()
                if not sequence:
                    raise ValueError(f"Empty FASTA record: {seq_id}")
                records.append((seq_id, sequence))

            header = line[1:].strip()
            seq_id = header.split()[0] if header else f"seq{len(records) + 1}"
            seq_lines = []
            continue

        if seq_id is None:
            raise ValueError("Invalid FASTA: sequence data found before header")
        seq_lines.append(line)

    if seq_id is not None:
        sequence = "".join(seq_lines).upper()
        if not sequence:
            raise ValueError(f"Empty FASTA record: {seq_id}")
        records.append((seq_id, sequence))

    if not records:
        raise ValueError("No FASTA records found")

    return records
