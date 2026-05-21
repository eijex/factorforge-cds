"""
FactorForge Command Line Interface
For Linux servers and automation
"""

import argparse
import json
import sys
from pathlib import Path

# Add src root to path
src_root = next(p for p in Path(__file__).resolve().parents if p.name == "src")
sys.path.insert(0, str(src_root))

def main():
    parser = argparse.ArgumentParser(description="FactorForge - AI-Powered Codon Optimization")

    parser.add_argument("input", help="Input DNA sequence or FASTA file")

    parser.add_argument("-o", "--output", help="Output file (JSON format)", default=None)

    parser.add_argument(
        "-m", "--model", help="Model checkpoint path", default="outputs/checkpoints/phase3"
    )

    parser.add_argument("--batch", action="store_true", help="Process as FASTA file (batch mode)")

    parser.add_argument("--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    try:
        from factorforge.engines.v1_archived.evaluation.metrics import BiologicalMetrics
        from factorforge.engines.v1_archived.tokenization.codon_tokenizer import CodonTokenizer
    except ImportError as exc:
        raise SystemExit(
            "v1 dependencies not installed. Install with: pip install -e \".[v1]\""
        ) from exc

    # Load tokenizer
    if args.verbose:
        print(f"Loading tokenizer from {args.model}...")

    tokenizer_path = Path(args.model) / "tokenizer"
    if not tokenizer_path.exists():
        print(f"Error: Tokenizer not found at {tokenizer_path}")
        return

    tokenizer = CodonTokenizer.load(str(tokenizer_path))

    # Get sequence
    sequences = []
    if Path(args.input).exists():
        # File input
        with open(args.input) as f:
            content = f.read()

        # Handle FASTA
        if content.startswith(">"):
            if args.batch:
                sequences = parse_fasta(content)
            else:
                sequences = [extract_first_sequence(content)]
        else:
            sequences = [content.strip()]
    else:
        # Direct sequence input
        sequences = [args.input]

    # Process sequences
    results = []
    bio = BiologicalMetrics()

    for i, seq in enumerate(sequences):
        # Clean sequence
        seq = "".join(c for c in seq.upper() if c in "ATGC")

        if len(seq) == 0:
            continue

        if len(seq) % 3 != 0:
            if args.verbose:
                print(f"Warning: Sequence {i+1} length not multiple of 3, skipping...")
            continue

        # Tokenize
        tokens = tokenizer.encode(seq)

        # Metrics
        metrics = bio.evaluate_sequence_quality(seq)
        is_quality, checks = bio.is_high_quality(seq)

        result = {
            "sequence_id": i + 1,
            "length": len(seq),
            "num_tokens": len(tokens),
            "compression_ratio": round(len(seq) / len(tokens), 2) if len(tokens) > 0 else 0,
            "gc_content": metrics["gc_content"],
            "cai": metrics["cai"],
            "rare_codon_freq": metrics["rare_codon_freq"],
            "high_quality": is_quality,
            "quality_checks": checks,
        }

        results.append(result)

        if args.verbose:
            print(f"\nSequence {i+1}:")
            print(f"  Length: {result['length']} bp")
            print(f"  Tokens: {result['num_tokens']}")
            print(f"  GC%: {result['gc_content']:.1f}%")
            print(f"  CAI: {result['cai']:.3f}")
            print(f"  Quality: {'✅ PASS' if is_quality else '❌ FAIL'}")

    # Output
    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\n✅ Results saved to {args.output}")
    else:
        print("\n" + json.dumps(results, indent=2))


def parse_fasta(content):
    """Parse FASTA file with multiple sequences"""
    sequences = []
    current_seq = []

    for line in content.split("\n"):
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq))
                current_seq = []
        else:
            current_seq.append(line.strip())

    if current_seq:
        sequences.append("".join(current_seq))

    return sequences


def extract_first_sequence(content):
    """Extract first sequence from FASTA"""
    lines = content.split("\n")
    sequence = []

    for line in lines:
        if line.startswith(">"):
            continue
        sequence.append(line.strip())

    return "".join(sequence)


if __name__ == "__main__":
    main()
