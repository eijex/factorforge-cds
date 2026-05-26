"""
generate_training_data.py — Phase 1: v2 엔진으로 학습 트리플렛 생성

입력:  data/raw/uniprot_plant_proteins.fasta  (Phase 0 출력)
출력:  data/training/training_triplets.jsonl

각 단백질 × 6 프로파일 = 트리플렛:
{protein_id, protein_seq, profile, dna_seq, cai, gc_percent, polya_signals}

Usage:
    PYTHONPATH=src python scripts/2_preprocessing/generate_training_data.py
    PYTHONPATH=src python scripts/2_preprocessing/generate_training_data.py --input data/raw/uniprot_plant_proteins.fasta --output data/training/training_triplets.jsonl
"""

import argparse
import json
import sys
import time
from pathlib import Path

PROFILES = ["balanced", "high_cai", "gc_target", "assembly_friendly", "ramp", "viral_delivery"]

DEFAULT_INPUT = Path("data/raw/uniprot_plant_proteins.fasta")
DEFAULT_OUTPUT = Path("data/training/training_triplets.jsonl")


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """FASTA → [(protein_id, sequence), ...] 반환. UniProt 헤더 파싱."""
    entries = []
    current_id = None
    current_seq = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    entries.append((current_id, "".join(current_seq)))
                # UniProt 헤더: >sp|P12345|GENE_ORG Description
                # 또는: >tr|A0A0A0|... / >P12345 ...
                header = line[1:]
                parts = header.split("|")
                if len(parts) >= 2:
                    current_id = parts[1]  # UniProt accession
                else:
                    current_id = header.split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_id and current_seq:
        entries.append((current_id, "".join(current_seq)))

    return entries


def generate(input_path: Path, output_path: Path) -> None:
    # v2 엔진 로드 (PYTHONPATH=src 필요)
    try:
        from factorforge.engines.registry import EngineRegistry
        optimizer = EngineRegistry.get("v2")
    except ImportError:
        print("[ERROR] factorforge 모듈을 찾을 수 없습니다.", file=sys.stderr)
        print("실행 방법: PYTHONPATH=src python scripts/2_preprocessing/generate_training_data.py", file=sys.stderr)
        sys.exit(1)

    proteins = parse_fasta(input_path)
    total = len(proteins)
    print(f"단백질 {total}개 로드 완료 → {total * len(PROFILES)}개 트리플렛 생성 예정\n")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    failed = 0
    start = time.time()

    with open(output_path, "w", encoding="utf-8") as out:
        for i, (protein_id, protein_seq) in enumerate(proteins, 1):
            for profile in PROFILES:
                try:
                    result = optimizer.optimize(protein_seq, profile=profile)

                    polya_signals = len(
                        result.metadata.get("scan_results", {}).get("polya", [])
                    )

                    record = {
                        "protein_id": protein_id,
                        "protein_seq": protein_seq,
                        "profile": profile,
                        "dna_seq": result.sequence,
                        "cai": round(result.metrics.get("cai", 0.0), 4),
                        "gc_percent": round(result.metrics.get("gc_percent", 0.0), 2),
                        "polya_signals": polya_signals,
                    }
                    out.write(json.dumps(record, ensure_ascii=False) + "\n")
                    written += 1

                except Exception as e:
                    failed += 1
                    print(f"  [SKIP] {protein_id} / {profile}: {e}", file=sys.stderr)

            # 진행 상황 출력 (50개마다)
            if i % 50 == 0 or i == total:
                elapsed = time.time() - start
                rate = written / elapsed if elapsed > 0 else 0
                eta = (total - i) * len(PROFILES) / rate if rate > 0 else 0
                print(
                    f"진행: {i}/{total} 단백질 ({written}개 트리플렛, "
                    f"{failed}개 실패, {elapsed:.0f}s 경과, ETA {eta:.0f}s)"
                )

    elapsed = time.time() - start
    print(f"\n완료: {written}개 트리플렛 저장 → {output_path}")
    print(f"실패: {failed}개 | 소요 시간: {elapsed:.0f}초 ({elapsed/60:.1f}분)")
    print(f"파일 크기: {output_path.stat().st_size / 1024 / 1024:.1f} MB")


def main() -> None:
    parser = argparse.ArgumentParser(description="v2 엔진으로 학습 트리플렛 생성 (Phase 1)")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if not args.input.exists():
        print(f"[ERROR] 입력 파일 없음: {args.input}", file=sys.stderr)
        print("Phase 0 먼저 실행: python scripts/1_data_preparation/download_uniprot.py", file=sys.stderr)
        sys.exit(1)

    generate(args.input, args.output)


if __name__ == "__main__":
    main()
