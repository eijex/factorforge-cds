"""Extract per-token ESM2 embeddings for FactorForge v3."""

from __future__ import annotations

import argparse
from pathlib import Path

import torch


def _protein_id(record_id: str) -> str:
    parts = record_id.split("|")
    if len(parts) >= 2:
        return parts[1]
    return record_id.split()[0]


def extract_embeddings(
    fasta_path: str,
    output_dir: str,
    model_name: str = "esm2_t6_8M_UR50D",
    batch_size: int = 16,
    max_length: int = 512,
) -> None:
    """Save one ``.pt`` file per protein with shape ``(seq_len, 320)``."""
    import esm
    from Bio import SeqIO

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    print(f"Loading ESM2 model: {model_name}")
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    batch_converter = alphabet.get_batch_converter()

    records = list(SeqIO.parse(fasta_path, "fasta"))
    to_process: list[tuple[str, str]] = []
    for record in records:
        protein_id = _protein_id(record.id)
        out_file = output / f"{protein_id}.pt"
        if not out_file.exists():
            sequence = str(record.seq)[:max_length]
            to_process.append((protein_id, sequence))

    print(f"Total: {len(records)} | To process: {len(to_process)}")

    for batch_start in range(0, len(to_process), batch_size):
        batch = to_process[batch_start : batch_start + batch_size]
        _, _, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[6], return_contacts=False)

        token_repr = results["representations"][6]

        for batch_idx, (protein_id, sequence) in enumerate(batch):
            per_token = token_repr[batch_idx, 1 : len(sequence) + 1, :].cpu()
            torch.save(
                {"sequence": sequence, "embeddings": per_token},
                output / f"{protein_id}.pt",
            )

        done = min(batch_start + batch_size, len(to_process))
        print(f"Processed: {done}/{len(to_process)}", end="\r")

    print(f"\nDone. Embeddings saved to {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", default="data/raw/uniprot_nbenthamiana_extended.fasta")
    parser.add_argument("--output", default="data/embeddings/per_token")
    parser.add_argument("--model", default="esm2_t6_8M_UR50D")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--max-length", type=int, default=512)
    args = parser.parse_args()
    extract_embeddings(
        args.fasta,
        args.output,
        args.model,
        args.batch_size,
        args.max_length,
    )


if __name__ == "__main__":
    main()
