"""
Codon Table Builder for FactorForge profile engine.
Build blended codon usage tables from multiple data sources for optimized CAI calculation.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def build_golden_set(
    high_expression_path: str | Path,
    empirical_path: str | Path,
    blend_ratio: float = 0.7,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """
    Build a 'golden set' codon table blending high-expression and empirical data.

    The golden set uses codon frequencies biased toward highly expressed genes,
    providing more accurate CAI reference weights per Sharp & Li (1987).

    Args:
        high_expression_path: Path to high-expression reference frequencies JSON.
            Must contain a "codon_usage" dict mapping amino acids to codon frequencies.
        empirical_path: Path to empirical codon table JSON (e.g. RNA-seq expression-weighted frequencies).
            Must contain "codons" and "amino_acids" sections.
        blend_ratio: Weight for high-expression data (0.0-1.0). Default 0.7.
        output_path: Optional path to write the blended table.

    Returns:
        Blended codon table dict in the standard format (codons + amino_acids).
    """
    if not 0.0 <= blend_ratio <= 1.0:
        raise ValueError(f"blend_ratio must be between 0.0 and 1.0, got {blend_ratio}")

    with open(high_expression_path, "r", encoding="utf-8") as f:
        high_expr = json.load(f)

    with open(empirical_path, "r", encoding="utf-8") as f:
        empirical = json.load(f)

    # Extract codon_usage from high-expression source
    high_usage = high_expr.get("codon_usage", {})

    # Extract per-codon frequencies from empirical source
    empirical_codons = empirical.get("codons", {})

    # Build blended frequencies per amino acid
    blended_codons: dict[str, dict[str, Any]] = {}
    blended_amino_acids: dict[str, dict[str, Any]] = {}

    for aa, high_codons in high_usage.items():
        aa_freqs: dict[str, float] = {}

        for codon, high_freq in high_codons.items():
            # Get empirical frequency
            emp_info = empirical_codons.get(codon, {})
            emp_freq = emp_info.get("frequency", high_freq) if emp_info else high_freq

            # Blend: weighted average
            blended_freq = blend_ratio * high_freq + (1 - blend_ratio) * emp_freq
            aa_freqs[codon] = blended_freq

        # Normalize per amino acid (frequencies must sum to 1.0)
        total = sum(aa_freqs.values())
        if total > 0:
            aa_freqs = {c: round(f / total, 4) for c, f in aa_freqs.items()}

        # Build codons section entries
        for codon, freq in aa_freqs.items():
            blended_codons[codon] = {
                "aa": aa,
                "frequency": freq,
                "per_thousand": round(freq * 1000 / len(aa_freqs), 1),
            }

        # Build amino_acids section
        sorted_codons = sorted(aa_freqs.keys(), key=lambda c: aa_freqs[c], reverse=True)
        blended_amino_acids[aa] = {
            "codons": sorted_codons,
            "preferred": sorted_codons[0] if sorted_codons else "",
        }

    result: dict[str, Any] = {
        "organism": high_expr.get("species", empirical.get("organism", "Unknown")),
        "source": (
            f"Golden Set ({int(blend_ratio * 100)}% high-expression "
            f"+ {int((1 - blend_ratio) * 100)}% empirical)"
        ),
        "blend_ratio": blend_ratio,
        "codons": blended_codons,
        "amino_acids": blended_amino_acids,
        "gc_content": empirical.get("gc_content", {"overall": 0.44}),
    }

    if output_path is not None:
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=4, ensure_ascii=False)
        logger.info(f"Golden set written to {output_path}")

    return result
