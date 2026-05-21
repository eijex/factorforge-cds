import sys

if sys.stdout.encoding != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except:
        pass

import json
import os
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.evaluation.metrics import Evaluator
from src.models.transformer import CodonTransformer, generate_square_subsequent_mask
from src.tokenization.codon_tokenizer import CodonTokenizer
from src.training.dataset import CodonDataset


def evaluate():
    print("=" * 70)
    print("🧬 Phase 4: Model Evaluation")
    print("=" * 70)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"💻 Device: {device}")

    # 1. Load Checkpoint
    checkpoint_path = Path("outputs/checkpoints/phase3/model.pt")
    tokenizer_path = Path("outputs/checkpoints/phase3/tokenizer")

    if not checkpoint_path.exists():
        print(f"❌ Checkpoint not found at {checkpoint_path}")
        return

    print(f"\n📂 Loading checkpoint from {checkpoint_path}...")
    checkpoint = torch.load(checkpoint_path, map_location=device)

    # 2. Load Tokenizer
    print("📚 Loading tokenizer...")
    tokenizer = CodonTokenizer.load(str(tokenizer_path))
    vocab_size = tokenizer.tokenizer.get_vocab_size()

    # 3. Load Benchmark Dataset
    benchmark_path = Path("data/benchmark/dataset.json")
    print(f"📊 Loading benchmark data from {benchmark_path}...")
    with open(benchmark_path, "r") as f:
        benchmark_data = json.load(f)

    sequences = [item["sequence"] for item in benchmark_data]
    labels = [item["label"] for item in benchmark_data]
    ids = [item["id"] for item in benchmark_data]

    dataset = CodonDataset(sequences, tokenizer, max_len=512)
    loader = DataLoader(dataset, batch_size=1, shuffle=False)

    # 4. Initialize Model
    print("🤖 Initializing model...")
    state_dict = checkpoint["model_state_dict"]

    # Recreate model with same config as training
    # Loading config from checkpoint would be better, but we hardcode for Phase 4 minimal
    model = CodonTransformer(
        vocab_size=vocab_size,
        d_model=128,
        nhead=4,
        num_layers=2,
        dim_feedforward=512,
        dropout=0.1,
        max_len=512,
        model_type="decoder",
    ).to(device)

    model.load_state_dict(state_dict)
    model.eval()

    # 5. Run Evaluation
    print("\n🚀 Running evaluation...")

    criterion = nn.CrossEntropyLoss(ignore_index=0, reduction="mean")
    results = []

    total_loss = 0
    total_acc = 0

    with torch.no_grad():
        for i, batch in enumerate(loader):
            input_tensor = batch.to(device)
            src = input_tensor[:, :-1]
            tgt = input_tensor[:, 1:]

            seq_len = src.size(1)
            mask = generate_square_subsequent_mask(seq_len).to(device)

            output = model(src, mask=mask)

            # Reshape
            output_flat = output.reshape(-1, vocab_size)
            tgt_flat = tgt.reshape(-1)

            # Metrics
            loss = criterion(output_flat, tgt_flat).item()
            ppl = Evaluator.calculate_perplexity(loss)
            acc = Evaluator.calculate_accuracy(output, tgt)
            gc = Evaluator.calculate_gc_content(sequences[i])

            total_loss += loss
            total_acc += acc

            results.append(
                {
                    "id": ids[i],
                    "label": labels[i],
                    "loss": round(loss, 4),
                    "perplexity": round(ppl, 4),
                    "accuracy": round(acc, 4),
                    "gc_content": round(gc, 2),
                }
            )

    avg_loss = total_loss / len(loader)
    avg_ppl = Evaluator.calculate_perplexity(avg_loss)
    avg_acc = total_acc / len(loader)

    print("\n📈 Aggregate Results:")
    print(f"  Average Loss: {avg_loss:.4f}")
    print(f"  Average Perplexity: {avg_ppl:.4f}")
    print(f"  Average Accuracy: {avg_acc:.4f}")

    # Segmented Analysis
    print("\n🔍 Analysis by Label:")
    positives = [r for r in results if r["label"] == "positive"]
    negatives = [r for r in results if r["label"] == "negative"]

    if positives:
        pos_ppl = sum(r["perplexity"] for r in positives) / len(positives)
        print(f"  Positive Samples PPL: {pos_ppl:.4f}")

    if negatives:
        neg_ppl = sum(r["perplexity"] for r in negatives) / len(negatives)
        print(f"  Negative Samples PPL: {neg_ppl:.4f}")

    # Save Results
    output_file = Path("outputs/evaluation/results.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    final_output = {
        "summary": {"avg_loss": avg_loss, "avg_perplexity": avg_ppl, "avg_accuracy": avg_acc},
        "details": results,
    }

    with open(output_file, "w") as f:
        json.dump(final_output, f, indent=2)

    print(f"\n✅ Results saved to {output_file}")


if __name__ == "__main__":
    evaluate()
