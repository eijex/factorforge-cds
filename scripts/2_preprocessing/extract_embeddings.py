"""
extract_embeddings.py — Phase 2: ESM2로 단백질 임베딩 추출

입력:  data/training/training_triplets.jsonl  (Phase 1 출력)
출력:  data/embeddings/protein_embeddings.pt

모델: facebook/esm2_t6_8M_UR50D  (소형, CPU 가능)
임베딩: mean pooling → (N, 320) 텐서
중복 단백질은 한 번만 계산 (1,000개 고유 서열)

Usage:
    pip install torch transformers
    python scripts/2_preprocessing/extract_embeddings.py

    # GPU 사용 시
    python scripts/2_preprocessing/extract_embeddings.py --device cuda

    # 더 큰 모델 (정확도↑, 메모리↑)
    python scripts/2_preprocessing/extract_embeddings.py --model facebook/esm2_t33_650M_UR50D
"""

import argparse
import json
import sys
import time
from pathlib import Path

DEFAULT_INPUT = Path("data/training/training_triplets.jsonl")
DEFAULT_OUTPUT = Path("data/embeddings/protein_embeddings.pt")
DEFAULT_MODEL = "facebook/esm2_t6_8M_UR50D"
BATCH_SIZE = 16


def load_unique_proteins(jsonl_path: Path) -> dict[str, str]:
    """JSONL에서 고유 단백질 추출 → {protein_id: protein_seq}"""
    proteins = {}
    with open(jsonl_path, encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)
            pid = record["protein_id"]
            if pid not in proteins:
                proteins[pid] = record["protein_seq"]
    return proteins


def extract(input_path: Path, output_path: Path, model_name: str, device: str) -> None:
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError:
        print("[ERROR] torch 또는 transformers 미설치", file=sys.stderr)
        print("설치: pip install torch transformers", file=sys.stderr)
        sys.exit(1)

    proteins = load_unique_proteins(input_path)
    protein_ids = list(proteins.keys())
    seqs = list(proteins.values())
    n = len(protein_ids)
    print(f"고유 단백질: {n}개")
    print(f"모델: {model_name}")
    print(f"디바이스: {device}\n")

    print("모델 로딩 중...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model = model.to(device)
    model.eval()

    import torch
    all_embeddings = []
    start = time.time()

    with torch.no_grad():
        for batch_start in range(0, n, BATCH_SIZE):
            batch_ids = protein_ids[batch_start:batch_start + BATCH_SIZE]
            batch_seqs = seqs[batch_start:batch_start + BATCH_SIZE]

            # 긴 서열 512aa로 자르기 (ESM2 최대 입력 길이)
            batch_seqs_trimmed = [s[:512] for s in batch_seqs]

            inputs = tokenizer(
                batch_seqs_trimmed,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=514,  # 512aa + 2 special tokens
            )
            inputs = {k: v.to(device) for k, v in inputs.items()}

            outputs = model(**inputs)

            # Mean pooling (패딩 제외)
            attention_mask = inputs["attention_mask"]
            token_embeddings = outputs.last_hidden_state  # (B, L, D)
            mask_expanded = attention_mask.unsqueeze(-1).float()
            sum_embeddings = (token_embeddings * mask_expanded).sum(dim=1)
            sum_mask = mask_expanded.sum(dim=1).clamp(min=1e-9)
            embeddings = sum_embeddings / sum_mask  # (B, D)

            all_embeddings.append(embeddings.cpu())

            done = min(batch_start + BATCH_SIZE, n)
            elapsed = time.time() - start
            print(f"진행: {done}/{n} ({elapsed:.0f}s)")

    all_embeddings = torch.cat(all_embeddings, dim=0)  # (N, D)
    print(f"\n임베딩 shape: {all_embeddings.shape}")

    # 저장: {protein_id → index} 매핑 + 텐서
    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "embeddings": all_embeddings,       # (N, D) float32
            "protein_ids": protein_ids,          # [str, ...]
            "model": model_name,
            "dim": all_embeddings.shape[1],
        },
        output_path,
    )

    elapsed = time.time() - start
    print(f"저장 완료 → {output_path}")
    print(f"파일 크기: {output_path.stat().st_size / 1024:.1f} KB | 소요: {elapsed:.0f}초")


def main() -> None:
    parser = argparse.ArgumentParser(description="ESM2 단백질 임베딩 추출 (Phase 2)")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--model", default=DEFAULT_MODEL, help=f"HuggingFace 모델 ID (기본: {DEFAULT_MODEL})")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda", "mps"], help="연산 디바이스 (기본: cpu)")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"[ERROR] 입력 파일 없음: {args.input}", file=sys.stderr)
        print("Phase 1 먼저 실행: python scripts/2_preprocessing/generate_training_data.py", file=sys.stderr)
        sys.exit(1)

    extract(args.input, args.output, args.model, args.device)


if __name__ == "__main__":
    main()
