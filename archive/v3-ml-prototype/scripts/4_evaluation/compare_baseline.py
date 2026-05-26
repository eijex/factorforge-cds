"""Compare FactorForge v2, v3 Run 1, and v3 Run 2 metrics."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "3_training"))

from factorforge.engines.v3.metrics import compute_cai, compute_gc, load_codon_usage_table


TEST_PROTEINS = {
    "GFP": (
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSR"
        "YPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYN"
        "SHVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVL"
        "LEFVTAAGITHGMDELYK"
    ),
    "CD47": (
        "MWPLVAALLLGSACCGSAQLLFNKTKSVEHSDGDLVNEVDGSNFTVSLEPGGRRITMQLKPKDGEFIQSPTR"
        "TLDQFTFVQLNESKEVEGMAYRMV"
    ),
    "mCherry": (
        "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMY"
        "GSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKK"
        "TMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIV"
        "EQYERAEGRHSTGGMDELYK"
    ),
}

V2_BASELINE = {
    "GFP": {"cai": 0.9093, "gc": 56.4, "enc": 42.22},
    "CD47": {"cai": 0.9010, "gc": 59.7, "enc": 42.99},
    "mCherry": {"cai": 0.9340, "gc": 58.2, "enc": 41.35},
}

V3_RUN1_BASELINE = {
    "GFP": {"cai": 0.9737, "gc": 60.9, "enc": None},
}

CODON_TO_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def format_metric(value: float | None, digits: int = 4, suffix: str = "") -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}{suffix}"


def family_homozygosity(counts: list[int]) -> float | None:
    total = sum(counts)
    if total <= 1:
        return None
    return ((sum(count * count for count in counts) - total) / (total * (total - 1)))


def compute_enc(dna_sequence: str) -> float | None:
    seq = dna_sequence.upper().replace("U", "T")
    aa_to_codons: dict[str, dict[str, int]] = {
        aa: {codon: 0 for codon in codons}
        for aa, codons in AA_TO_CODONS.items()
        if len(codons) in {2, 3, 4, 6}
    }
    for index in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[index : index + 3]
        aa = CODON_TO_AA.get(codon)
        if aa is None:
            continue
        if aa in aa_to_codons:
            aa_to_codons[aa][codon] += 1

    by_degeneracy: dict[int, list[float]] = {2: [], 3: [], 4: [], 6: []}
    for codon_counts in aa_to_codons.values():
        degeneracy = len(codon_counts)
        if degeneracy not in by_degeneracy:
            continue
        homozygosity = family_homozygosity(list(codon_counts.values()))
        if homozygosity and homozygosity > 0:
            by_degeneracy[degeneracy].append(homozygosity)

    means = {
        degeneracy: (sum(values) / len(values) if values else None)
        for degeneracy, values in by_degeneracy.items()
    }
    if any(value is None or value <= 0 for value in means.values()):
        return None
    enc = 2 + 9 / means[2] + 1 / means[3] + 5 / means[4] + 3 / means[6]
    return min(61.0, max(20.0, float(enc))) if math.isfinite(enc) else None


def resolve_checkpoint(path: str | None) -> Path | None:
    if not path:
        return None
    checkpoint = Path(path)
    if checkpoint.is_dir():
        model_file = checkpoint / "pytorch_model.pt"
        if model_file.exists():
            return model_file
        step_files = sorted(checkpoint.glob("step_*.pt"))
        if step_files:
            return step_files[-1]
    return checkpoint if checkpoint.exists() else None


def generate_run2_sequences(
    checkpoint_path: Path,
    config_path: Path,
    protein_names: list[str],
) -> dict[str, str]:
    import torch
    import yaml
    from esm import pretrained
    from train_v3_esm2_bart import build_model

    from factorforge.engines.v3.tokenizer import CodonTokenizer

    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    codon_tok = CodonTokenizer.default()
    model = build_model(cfg).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=True)
    state_dict = checkpoint.get("model", checkpoint) if isinstance(checkpoint, dict) else checkpoint
    model.load_state_dict(state_dict)
    model.eval()

    esm_model, alphabet = pretrained.load_model_and_alphabet(cfg["esm2"]["model_name"])
    esm_model = esm_model.to(device)
    esm_model.eval()
    batch_converter = alphabet.get_batch_converter()

    generated: dict[str, str] = {}
    for name in protein_names:
        sequence = TEST_PROTEINS[name]
        _, _, tokens = batch_converter([(name, sequence)])
        tokens = tokens.to(device)
        with torch.no_grad():
            embeddings = esm_model(tokens, repr_layers=[6], return_contacts=False)["representations"][6]
            per_token = embeddings[:, 1 : len(sequence) + 1, :]
            output_ids = model.generate(
                per_token,
                max_new_tokens=len(sequence) + 1,
                bos_token_id=codon_tok.bos_token_id,
                eos_token_id=codon_tok.eos_token_id,
                pad_token_id=codon_tok.pad_token_id,
            )
        generated[name] = codon_tok.decode(output_ids[0].tolist())[: len(sequence) * 3]
    return generated


def row_for_metrics(protein: str, method: str, metrics: dict[str, float | None]) -> list[str]:
    return [
        protein,
        method,
        format_metric(metrics.get("cai"), digits=4),
        format_metric(metrics.get("gc"), digits=1, suffix="%"),
        format_metric(metrics.get("enc"), digits=2),
    ]


def print_table(rows: list[list[str]]) -> None:
    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    for row_index, row in enumerate(rows):
        print(" | ".join(value.ljust(widths[index]) for index, value in enumerate(row)))
        if row_index == 0:
            print("-|-".join("-" * width for width in widths))


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare v2, v3 Run 1, and v3 Run 2 baselines.")
    parser.add_argument("checkpoint_path", nargs="?", help="Run 2 checkpoint file or directory.")
    parser.add_argument("--checkpoint", dest="checkpoint_option", help="Run 2 checkpoint file or directory.")
    parser.add_argument("--config", default="configs/v3_training_config_run2.yml")
    parser.add_argument("--proteins", nargs="+", default=["GFP", "CD47", "mCherry"])
    parser.add_argument("--baseline", default=None, help="Accepted for job-run compatibility.")
    args = parser.parse_args()

    proteins = [name for name in args.proteins if name in TEST_PROTEINS]
    if not proteins:
        raise ValueError("No supported proteins selected. Use GFP, CD47, and/or mCherry.")

    checkpoint = resolve_checkpoint(args.checkpoint_option or args.checkpoint_path)
    run2_sequences: dict[str, str] = {}
    if checkpoint:
        run2_sequences = generate_run2_sequences(checkpoint, Path(args.config), proteins)

    rows = [["Protein", "Method", "CAI", "GC%", "ENC"]]
    table = load_codon_usage_table(ROOT / "data" / "nbenthamiana_codons.json")
    for protein in proteins:
        rows.append(row_for_metrics(protein, "v2", V2_BASELINE[protein]))
        rows.append(row_for_metrics(protein, "v3 Run 1", V3_RUN1_BASELINE.get(protein, {})))
        if protein in run2_sequences:
            dna = run2_sequences[protein]
            rows.append(
                row_for_metrics(
                    protein,
                    "v3 Run 2",
                    {
                        "cai": compute_cai(dna, table),
                        "gc": compute_gc(dna),
                        "enc": compute_enc(dna),
                    },
                )
            )
        else:
            rows.append(row_for_metrics(protein, "v3 Run 2", {}))

    print_table(rows)
    if not checkpoint:
        print("\nRun 2 checkpoint not found/provided; v3 Run 2 metrics left blank.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
