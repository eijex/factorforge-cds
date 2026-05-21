"""
download_uniprot.py — Phase 0: UniProt 식물 단백질 서열 다운로드

Usage:
    python scripts/1_data_preparation/download_uniprot.py
    python scripts/1_data_preparation/download_uniprot.py --count 500 --output data/raw/custom.fasta

필터:
    - 분류: Viridiplantae (taxonomy_id:33090)
    - 길이: 50~500 aa
    - 데이터셋: Swiss-Prot (reviewed:true)

출력:
    data/raw/uniprot_plant_proteins.fasta
"""

import argparse
import sys
import time
from pathlib import Path

import requests

# UniProt REST API
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"

# Viridiplantae (taxonomy:33090), Swiss-Prot reviewed, 50~500aa
QUERY = "(taxonomy_id:33090) AND (reviewed:true) AND (length:[50 TO 500])"

DEFAULT_OUTPUT = Path("data/raw/uniprot_plant_proteins.fasta")
DEFAULT_COUNT = 1000
PAGE_SIZE = 500  # UniProt 최대 페이지 크기


def fetch_page(session: requests.Session, cursor: str | None, size: int) -> tuple[list[str], str | None]:
    """UniProt 한 페이지 다운로드. (entries, next_cursor) 반환."""
    params = {
        "query": QUERY,
        "format": "fasta",
        "size": size,
    }
    if cursor:
        params["cursor"] = cursor

    resp = session.get(UNIPROT_API, params=params, timeout=60)
    resp.raise_for_status()

    # Link 헤더에서 next cursor 추출
    next_cursor = None
    link_header = resp.headers.get("Link", "")
    if 'rel="next"' in link_header:
        # 형식: <https://...?cursor=XYZ&...>; rel="next"
        for part in link_header.split(","):
            if 'rel="next"' in part:
                url_part = part.strip().split(";")[0].strip("<> ")
                for token in url_part.split("&"):
                    if token.startswith("cursor="):
                        next_cursor = token[len("cursor="):]
                        break

    # FASTA 엔트리 분리
    text = resp.text
    entries = []
    current = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current:
                entries.append("\n".join(current))
            current = [line]
        elif line.strip():
            current.append(line)
    if current:
        entries.append("\n".join(current))

    return entries, next_cursor


def download(count: int, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "FactorForge/2.5.2 (eijex; research)"})

    collected: list[str] = []
    cursor: str | None = None
    page = 0

    print(f"UniProt 다운로드 시작: {count}개 목표")
    print(f"쿼리: {QUERY}\n")

    while len(collected) < count:
        remaining = count - len(collected)
        page_size = min(PAGE_SIZE, remaining)
        page += 1

        print(f"페이지 {page} 요청 중... (현재 {len(collected)}개 수집)")

        try:
            entries, cursor = fetch_page(session, cursor, page_size)
        except requests.HTTPError as e:
            print(f"[ERROR] HTTP 오류: {e}", file=sys.stderr)
            break
        except requests.RequestException as e:
            print(f"[ERROR] 네트워크 오류: {e}", file=sys.stderr)
            break

        if not entries:
            print("더 이상 결과 없음, 다운로드 완료.")
            break

        collected.extend(entries)
        print(f"  → {len(entries)}개 수신 (누계: {len(collected)}개)")

        if cursor is None:
            print("마지막 페이지 도달.")
            break

        # Rate limiting: UniProt 권장 (초당 1~2 요청)
        time.sleep(0.5)

    # 목표 수 초과분 제거
    collected = collected[:count]

    # FASTA 파일 저장
    with open(output, "w", encoding="utf-8") as f:
        f.write("\n\n".join(collected))
        f.write("\n")

    print(f"\n완료: {len(collected)}개 서열 저장 → {output}")
    print(f"파일 크기: {output.stat().st_size / 1024:.1f} KB")


def main() -> None:
    parser = argparse.ArgumentParser(description="UniProt 식물 단백질 FASTA 다운로드 (Phase 0)")
    parser.add_argument("--count", type=int, default=DEFAULT_COUNT, help=f"다운로드할 서열 수 (기본: {DEFAULT_COUNT})")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"출력 경로 (기본: {DEFAULT_OUTPUT})")
    args = parser.parse_args()

    download(args.count, args.output)


if __name__ == "__main__":
    main()
