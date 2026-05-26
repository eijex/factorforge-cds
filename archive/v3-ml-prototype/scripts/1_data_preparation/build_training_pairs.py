"""
Build v3 training pairs: protein sequence → v2-optimized codon sequence.

For each protein in the FASTA file:
  1. Run FactorForge v2 (balanced profile) to get the target codon sequence
  2. Check that a per-token ESM2 embedding exists
  3. Save {protein_id, sequence, codon_sequence} to training_pairs_v3.jsonl

Output format (one JSON per line):
    {"protein_id": "sp|P12345|...", "sequence": "MVSK...", "codon_sequence": "ATG GTC AGC ..."}

The codon_sequence uses space-separated triplets for tokenizer compatibility.

Usage:
    python scripts/1_data_preparation/build_training_pairs.py \\
        --fasta data/raw/uniprot_nbenthamiana_extended.fasta \\
        --embeddings data/embeddings/per_token \\
        --output data/training/training_pairs_v3.jsonl \\
        --profile balanced

Note:
    Run extract_per_token_esm2.py first to populate data/embeddings/per_token/.
    If --embeddings-required is set, only proteins with existing embeddings are included.
    Otherwise all proteins are included (embeddings extracted on-the-fly during training).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from Bio import SeqIO

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from factorforge.engines.v2.pipeline import OptimizationPipeline


STOP_CODONS = {"TAA", "TAG", "TGA"}


def dna_to_spaced_codons(dna: str) -> str:
    """Convert a DNA string to space-separated codon triplets."""
    return " ".join(dna[i : i + 3] for i in range(0, len(dna) - 2, 3))


def build_pairs(
    fasta_path: str,
    embeddings_dir: str,
    output_path: str,
    profile: str = "balanced",
    embeddings_required: bool = False,
    max_length: int = 512,
    skip_on_error: bool = True,
) -> None:
    fasta = Path(fasta_path)
    emb_dir = Path(embeddings_dir)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    pipeline = OptimizationPipeline(profile=profile)

    records = list(SeqIO.parse(fasta, "fasta"))
    print(f"Proteins in FASTA: {len(records)}")

    written = 0
    skipped_long = 0
    skipped_error = 0
    skipped_no_emb = 0

    with open(out, "w", encoding="utf-8") as f:
        for rec in records:
            protein_id = rec.id
            sequence = str(rec.seq).upper().replace("*", "")

            # Truncate to max_length
            if len(sequence) > max_length:
                sequence = sequence[:max_length]
                skipped_long += 1  # still included, just truncated

            # Check embedding exists
            emb_path = emb_dir / f"{protein_id}.pt"
            if embeddings_required and not emb_path.exists():
                skipped_no_emb += 1
                continue

            # Generate codon sequence via v2
            try:
                result = pipeline.run(sequence)
                dna = result.sequence

                # Validate: must be 3x protein length
                if len(dna) != len(sequence) * 3:
                    raise ValueError(f"Length mismatch: {len(dna)} != {len(sequence) * 3}")

                # Convert to space-separated codons (exclude stop codon)
                codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]
                codon_sequence = " ".join(c for c in codons if c not in STOP_CODONS)

                row = {
                    "protein_id": protein_id,
                    "sequence": sequence,
                    "codon_sequence": codon_sequence,
                }
                f.write(json.dumps(row, ensure_ascii=False) + "\n")
                written += 1

            except Exception as e:
                if skip_on_error:
                    skipped_error += 1
                else:
                    raise

            if written % 500 == 0 and written > 0:
                print(f"  Written: {written}", end="\r")

    print(f"\nDone.")
    print(f"  Written:          {written}")
    print(f"  Truncated (>={max_length} AA): {skipped_long}")
    print(f"  Skipped (error):  {skipped_error}")
    print(f"  Skipped (no emb): {skipped_no_emb}")
    print(f"  Output: {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build v3 training pairs from FASTA.")
    parser.add_argument(
        "--fasta",
        default="data/raw/uniprot_nbenthamiana_extended.fasta",
        help="Input FASTA file",
    )
    parser.add_argument(
        "--embeddings",
        default="data/embeddings/per_token",
        help="Directory with per-token ESM2 embeddings (.pt files)",
    )
    parser.add_argument(
        "--output",
        default="data/training/training_pairs_v3.jsonl",
        help="Output JSONL file",
    )
    parser.add_argument(
        "--profile",
        default="balanced",
        choices=["balanced", "high_cai", "gc_target", "assembly_friendly", "ramp", "viral_delivery"],
        help="FactorForge v2 optimization profile",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=512,
        help="Truncate proteins longer than this",
    )
    parser.add_argument(
        "--embeddings-required",
        action="store_true",
        help="Skip proteins without pre-computed ESM2 embeddings",
    )
    args = parser.parse_args()

    build_pairs(
        fasta_path=args.fasta,
        embeddings_dir=args.embeddings,
        output_path=args.output,
        profile=args.profile,
        max_length=args.max_length,
        embeddings_required=args.embeddings_required,
    )
