"""Small native FASTA reader/writer with privacy-safe output headers."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

from .validation import validate_sequence

HEADER_ALLOWLIST = ("engine", "host_profile", "profile", "sequence_hash")
BLOCKED_HEADER_TERMS = (
    "plantform",
    "confidential",
    "private",
    "secret",
    "partner",
    "yield",
    "wet-lab",
    "wet_lab",
    "clinical",
)
RAW_SEQUENCE_PATTERN = re.compile(r"[ACGTRYSWKMBDHVN]{20,}", re.IGNORECASE)
SAFE_HEADER_VALUE = re.compile(r"^[A-Za-z0-9_.:@/+ -]+$")


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    sequence: str
    metadata: Mapping[str, str] = field(default_factory=dict)


def _validate_header_value(value: str) -> str:
    normalized = str(value).strip()
    lowered = normalized.lower()
    if not normalized or len(normalized) > 120:
        raise ValueError("FASTA header values must contain 1-120 characters")
    if any(term in lowered for term in BLOCKED_HEADER_TERMS):
        raise ValueError("FASTA header contains blocked private or claim-related metadata")
    if RAW_SEQUENCE_PATTERN.search(normalized):
        raise ValueError("FASTA header must not contain a raw sequence")
    if not SAFE_HEADER_VALUE.fullmatch(normalized):
        raise ValueError("FASTA header contains unsupported characters")
    return normalized


def build_fasta_header(identifier: str, metadata: Mapping[str, object] | None = None) -> str:
    """Build an allowlist-only public FASTA header."""
    parts = [_validate_header_value(identifier)]
    for key in HEADER_ALLOWLIST:
        if metadata is None or key not in metadata or metadata[key] is None:
            continue
        value = _validate_header_value(str(metadata[key]))
        parts.append(f"{key}={value}")
    return " ".join(parts)


def parse_fasta(text: str, validation_mode: str | None = None) -> list[FastaRecord]:
    """Parse FASTA text, optionally validating each sequence alphabet."""
    records: list[FastaRecord] = []
    identifier: str | None = None
    sequence_lines: list[str] = []

    def append_record() -> None:
        if identifier is None:
            return
        sequence = "".join(sequence_lines)
        if validation_mode is not None:
            sequence = validate_sequence(sequence, validation_mode)
        records.append(FastaRecord(identifier=identifier, sequence=sequence))

    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            append_record()
            identifier = line[1:].strip()
            if not identifier:
                raise ValueError(f"FASTA header at line {line_no} is empty")
            sequence_lines = []
        elif identifier is None:
            raise ValueError(f"FASTA sequence found before first header at line {line_no}")
        else:
            sequence_lines.append(line)

    append_record()
    if not records:
        raise ValueError("FASTA input contains no records")
    return records


def format_fasta(
    records: Iterable[FastaRecord],
    *,
    validation_mode: str = "dna_strict",
    line_width: int = 60,
) -> str:
    """Serialize records using privacy-safe headers and validated sequences."""
    if line_width < 1:
        raise ValueError("line_width must be positive")

    lines: list[str] = []
    for record in records:
        header = build_fasta_header(record.identifier, record.metadata)
        sequence = validate_sequence(record.sequence, validation_mode)
        lines.append(f">{header}")
        lines.extend(sequence[i : i + line_width] for i in range(0, len(sequence), line_width))
    if not lines:
        raise ValueError("At least one FASTA record is required")
    return "\n".join(lines) + "\n"


def read_fasta(path: str | Path, validation_mode: str | None = None) -> list[FastaRecord]:
    return parse_fasta(Path(path).read_text(encoding="utf-8"), validation_mode)


def write_fasta(
    path: str | Path,
    records: Iterable[FastaRecord],
    *,
    validation_mode: str = "dna_strict",
    line_width: int = 60,
) -> Path:
    output_path = Path(path)
    output_path.write_text(
        format_fasta(records, validation_mode=validation_mode, line_width=line_width),
        encoding="utf-8",
    )
    return output_path
