# factorforge/src/factorforge/profiles/profile_snapshot.py
"""Immutable Host Profile Snapshot Loader and Validator."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class HostProfileSnapshot:
    """Immutable, content-addressed host codon and biological profile."""
    organism: str
    snapshot_version: str
    description: str
    codon_table: Dict[str, float]
    frequency_table: Dict[str, float]
    metadata: Dict[str, Any]
    digest: str

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> HostProfileSnapshot:
        """Construct profile snapshot from dictionary and verify digest."""
        organism = data["organism"]
        snapshot_version = data["snapshot_version"]
        description = data.get("description", "")
        codon_table = data["codon_table"]
        frequency_table = data.get("frequency_table", {})
        metadata = data.get("metadata", {})
        expected_digest = data.get("digest")

        calculated_digest = cls.compute_profile_digest(
            organism=organism,
            snapshot_version=snapshot_version,
            codon_table=codon_table,
            frequency_table=frequency_table,
        )

        if expected_digest and expected_digest != calculated_digest:
            raise ValueError(
                f"HostProfileSnapshot integrity mismatch! Expected {expected_digest}, computed {calculated_digest}"
            )

        return cls(
            organism=organism,
            snapshot_version=snapshot_version,
            description=description,
            codon_table=codon_table,
            frequency_table=frequency_table,
            metadata=metadata,
            digest=calculated_digest,
        )

    @classmethod
    def compute_profile_digest(
        cls,
        organism: str,
        snapshot_version: str,
        codon_table: Dict[str, float],
        frequency_table: Dict[str, float],
    ) -> str:
        """Compute deterministic SHA-256 fingerprint for host profile payload."""
        payload = {
            "organism": organism,
            "snapshot_version": snapshot_version,
            "codon_table": dict(sorted(codon_table.items())),
            "frequency_table": dict(sorted(frequency_table.items())),
        }
        canonical_json = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return f"sha256:{hashlib.sha256(canonical_json).hexdigest()}"

    def export_canonical_json(self) -> str:
        """Export snapshot as canonical JSON string."""
        payload = {
            "organism": self.organism,
            "snapshot_version": self.snapshot_version,
            "description": self.description,
            "codon_table": dict(sorted(self.codon_table.items())),
            "frequency_table": dict(sorted(self.frequency_table.items())),
            "metadata": self.metadata,
            "digest": self.digest,
        }
        return json.dumps(payload, indent=2, sort_keys=True)
