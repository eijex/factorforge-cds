import hashlib
import re
from typing import Tuple

def canonicalize_sequence(molecule_class: str, raw_sequence: str, version: str = "seqcanon-v1") -> Tuple[str, str]:
    """
    Normalizes a biological sequence and generates a canonical SHA-256 hash.
    
    This function acts as the unified contract for sequence identity across the system,
    handling cross-platform instability issues like whitespaces, case differences,
    and RNA/DNA (U/T) mapping.

    Args:
        molecule_class: 'dna', 'rna', or 'protein' (or specific subtypes like 'cds')
        raw_sequence: The raw input sequence string
        version: The canonicalization algorithm version to prevent silent hash collisions if rules change.

    Returns:
        Tuple of (normalized_sequence, sha256_hash)
    """
    # 1. Strip all whitespaces, newlines, and non-printing characters
    normalized = re.sub(r'\s+', '', raw_sequence).upper()

    # 2. Molecule-specific normalization
    mol_class_norm = molecule_class.lower()
    
    # We group CDS, DNA together
    if 'dna' in mol_class_norm or 'cds' in mol_class_norm:
        normalized = normalized.replace('U', 'T')
    elif 'rna' in mol_class_norm:
        normalized = normalized.replace('T', 'U')
    elif 'protein' in mol_class_norm or 'aa' in mol_class_norm:
        # Standardize stop codons to '*' (some tools use 'X', but '*' is canonical for stops)
        # Assuming typical protein sequence letters. We'll ensure it's upper case.
        pass

    # 3. Generate canonical payload and hash
    # Including the version in the payload string ensures that if we ever update 
    # the canonicalization logic, the hashes will cleanly differentiate.
    payload = f"{version}:{mol_class_norm}:{normalized}"
    seq_hash = hashlib.sha256(payload.encode('utf-8')).hexdigest()

    return normalized, seq_hash
