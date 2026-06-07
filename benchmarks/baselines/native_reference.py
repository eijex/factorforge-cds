"""Native reference anchor (NOT an optimizer). Returns native CDS as-is."""
from __future__ import annotations


def native_reference_cds(native_cds: str) -> str:
    """Identity: the benchmark's biological anchor is the native sequence itself."""
    return native_cds
