import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

if sys.stdout.encoding != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except:
        pass

from src.tokenization.baselines import NucleotideTokenizer, StandardCodonTokenizer
from src.tokenization.codon_tokenizer import CodonTokenizer


def main():
    print("=" * 70)
    print("⚖️  Phase 4: Tokenizer Comparison")
    print("=" * 70)

    # 1. Load Tokenizers
    print("📚 Loading tokenizers...")

    # FactorForge BPE
    bpe_path = Path("outputs/checkpoints/phase3/tokenizer")
    if not bpe_path.exists():
        print("❌ BPE Tokenizer not found! Please run training first.")
        return

    bpe_tokenizer = CodonTokenizer.load(str(bpe_path))

    # Baselines
    nuc_tokenizer = NucleotideTokenizer()
    std_tokenizer = StandardCodonTokenizer()

    print(f"  ✅ Nucleotide Vocab: {len(nuc_tokenizer.vocab)}")
    print(f"  ✅ Standard Codon Vocab: {len(std_tokenizer.vocab)}")
    print(f"  ✅ FactorForge BPE Vocab: {bpe_tokenizer.tokenizer.get_vocab_size()}")

    # 2. Load Benchmark Data
    data_path = Path("data/benchmark/dataset.json")
    with open(data_path, "r") as f:
        data = json.load(f)

    print(f"\n📊 Analyzing {len(data)} sequences from benchmark...")

    results = []

    for item in data:
        seq = item["sequence"]
        seq_len_bp = len(seq)

        # Encode
        enc_nuc = nuc_tokenizer.encode(seq)
        enc_std = std_tokenizer.encode(seq)
        enc_bpe = bpe_tokenizer.encode(seq)

        # Metrics
        len_nuc = len(enc_nuc)
        len_std = len(enc_std)
        len_bpe = len(enc_bpe)

        # Compression Ratio (BP / Token)
        # Higher is better (more info per token)
        ratio_nuc = seq_len_bp / len_nuc if len_nuc > 0 else 0
        ratio_std = seq_len_bp / len_std if len_std > 0 else 0
        ratio_bpe = seq_len_bp / len_bpe if len_bpe > 0 else 0

        results.append(
            {
                "id": item["id"],
                "type": item["type"],
                "len_bp": seq_len_bp,
                "len_nuc": len_nuc,
                "len_std": len_std,
                "len_bpe": len_bpe,
                "ratio_std": round(ratio_std, 2),
                "ratio_bpe": round(ratio_bpe, 2),
            }
        )

    # 3. Aggregate Analysis
    avg_len_nuc = np.mean([r["len_nuc"] for r in results])
    avg_len_std = np.mean([r["len_std"] for r in results])
    avg_len_bpe = np.mean([r["len_bpe"] for r in results])

    avg_ratio_std = np.mean([r["ratio_std"] for r in results])
    avg_ratio_bpe = np.mean([r["ratio_bpe"] for r in results])

    print("\n📈 Comparison Results (Averages):")
    print(f"  Original Length (BP): {np.mean([r['len_bp'] for r in results]):.1f}")
    print(f"  Nucleotide Tokens:    {avg_len_nuc:.1f} (Ratio: 1.0)")
    print(f"  Standard Codon Tokens:{avg_len_std:.1f} (Ratio: {avg_ratio_std:.2f})")
    print(f"  FactorForge BPE Tokens:{avg_len_bpe:.1f} (Ratio: {avg_ratio_bpe:.2f})")

    improvement = ((avg_len_std - avg_len_bpe) / avg_len_std) * 100
    print(f"\n🚀 BPE Compression Improvement over Standard: {improvement:.2f}%")

    # Save detailed results
    output_path = Path("outputs/evaluation/comparison_results.json")
    with open(output_path, "w") as f:
        json.dump(
            {
                "summary": {
                    "avg_ratio_std": avg_ratio_std,
                    "avg_ratio_bpe": avg_ratio_bpe,
                    "improvement_percent": improvement,
                },
                "details": results,
            },
            f,
            indent=2,
        )

    print(f"\n✅ Detailed results saved to {output_path}")


if __name__ == "__main__":
    main()
