"""Codon Distribution Policy definitions and validation."""

from __future__ import annotations

import hashlib
import json
import math

from pydantic import BaseModel, Field, model_validator

# Standard genetic code dictionary
CODON_TO_AA: dict[str, str] = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}

AA_TO_CODONS: dict[str, list[str]] = {}
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)


class CodonDistributionPolicy(BaseModel):
    """Specification of target synonymous codon distributions for specified residues."""

    schema_version: str = "0.1"
    policy_id: str
    host: str = "Nicotiana benthamiana"
    source_type: str = Field(
        default="expert_heuristic",
        description="Origin of distribution: 'expert_heuristic', 'custom', or 'learned'",
    )
    source_description: str = Field(
        default="",
        description="User-provided provenance annotation",
    )
    scope_region: str = Field(
        default="mature_chain",
        description="'mature_chain' or 'full_cds'",
    )
    preserve_signal_peptide: bool = True
    signal_peptide_aa_len: int = 20
    distributions: dict[str, dict[str, float]] = Field(
        ...,
        description="Target frequencies per amino acid, e.g. {'R': {'AGA': 0.5, 'AGG': 0.5}}",
    )
    allocation_method: str = "largest_remainder"
    tie_break: str = "lexical_codon_order"

    @model_validator(mode="after")
    def validate_distributions(self) -> CodonDistributionPolicy:
        if self.schema_version != "0.1":
            raise ValueError("Unsupported policy schema version")
        if self.scope_region not in {"mature_chain", "full_cds"}:
            raise ValueError("Unsupported policy region")
        if self.signal_peptide_aa_len < 0:
            raise ValueError("Signal peptide length cannot be negative")
        if self.allocation_method != "largest_remainder" or self.tie_break != "lexical_codon_order":
            raise ValueError("Unsupported allocation contract")
        if not self.distributions:
            raise ValueError("At least one distribution is required")
        for aa, codon_weights in self.distributions.items():
            aa_upper = aa.upper()
            if aa != aa_upper or aa == "*":
                raise ValueError("Use uppercase protein amino acid keys")
            valid_codons = AA_TO_CODONS.get(aa_upper, [])
            if not valid_codons:
                raise ValueError(f"Unknown or non-proteinogenic amino acid: '{aa}'")

            total_weight = 0.0
            for codon, weight in codon_weights.items():
                codon_upper = codon.upper()
                if codon_upper not in valid_codons:
                    raise ValueError(
                        f"Codon '{codon_upper}' does not code for amino acid '{aa_upper}' (valid: {valid_codons})"
                    )
                if codon != codon_upper or not math.isfinite(weight) or weight < 0.0:
                    raise ValueError(f"Codon weight cannot be negative: {codon_upper}={weight}")
                total_weight += weight

            if abs(total_weight - 1.0) > 1e-4:
                raise ValueError(
                    f"Codon weights for amino acid '{aa_upper}' must sum to 1.0 (got {total_weight:.6f})"
                )

        return self

    @property
    def canonical_sha256(self) -> str:
        """Deterministic canonical SHA-256 fingerprint of the policy configuration."""
        data = {
            "schema_version": self.schema_version,
            "policy_id": self.policy_id,
            "host": self.host,
            "source_type": self.source_type,
            "scope_region": self.scope_region,
            "preserve_signal_peptide": self.preserve_signal_peptide,
            "signal_peptide_aa_len": self.signal_peptide_aa_len,
            "distributions": {
                aa: {c: round(w, 6) for c, w in sorted(codons.items())}
                for aa, codons in sorted(self.distributions.items())
            },
            "allocation_method": self.allocation_method,
            "tie_break": self.tie_break,
        }
        encoded = json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()


# Standard presets according to Job 350 specification
DOUG_BALANCED_PRESET = CodonDistributionPolicy(
    policy_id="doug_balanced_reference_v1",
    host="Nicotiana benthamiana",
    source_type="expert_heuristic",
    source_description="Doug feedback initial reference profile (Arg 1:1, Ser 1:1:1:1:1)",
    scope_region="mature_chain",
    preserve_signal_peptide=True,
    signal_peptide_aa_len=20,
    distributions={
        "R": {"AGA": 0.50, "AGG": 0.50},
        "S": {"TCT": 0.20, "TCC": 0.20, "TCA": 0.20, "AGT": 0.20, "AGC": 0.20},
    },
)

ILLUSTRATIVE_SKEWED_PRESET = CodonDistributionPolicy(
    policy_id="illustrative_skewed_v1",
    host="Nicotiana benthamiana",
    source_type="custom",
    source_description="Illustrative skewed exploratory distribution",
    scope_region="mature_chain",
    preserve_signal_peptide=True,
    signal_peptide_aa_len=20,
    distributions={
        "R": {"AGA": 0.70, "AGG": 0.30},
        "S": {"TCT": 0.35, "TCC": 0.25, "TCA": 0.15, "AGT": 0.15, "AGC": 0.10},
    },
)

FACTORFORGE_RECOMMENDED_PRESET = CodonDistributionPolicy(
    policy_id="factorforge_host_optimized_v1",
    host="Nicotiana benthamiana",
    source_type="custom",
    source_description="Illustrative experimental distribution; not learned, host-calibrated, or validated",
    scope_region="mature_chain",
    preserve_signal_peptide=True,
    signal_peptide_aa_len=20,
    distributions={
        "R": {"AGA": 0.49, "AGG": 0.51},
        "S": {"TCT": 0.20, "TCC": 0.24, "TCA": 0.13, "AGT": 0.17, "AGC": 0.26},
    },
)
