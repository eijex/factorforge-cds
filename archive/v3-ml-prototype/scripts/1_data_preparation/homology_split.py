"""Split v3 Run 2 training pairs by homology clusters when CD-HIT is available."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


DEFAULT_INPUT = Path("data/training/training_pairs_v3_run2.jsonl")
DEFAULT_TRAIN = Path("data/training/train_v3_run2.jsonl")
DEFAULT_EVAL = Path("data/training/eval_v3_run2.jsonl")


def load_pairs(path: Path) -> list[dict[str, Any]]:
    pairs: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                pairs.append(json.loads(line))
    return pairs


def pair_id(pair: dict[str, Any]) -> str:
    value = pair.get("id") or pair.get("protein_id")
    if not isinstance(value, str) or not value:
        raise ValueError("Each pair must include a non-empty 'id' or 'protein_id'")
    return value


def pair_aa(pair: dict[str, Any]) -> str:
    value = pair.get("aa_seq") or pair.get("sequence")
    if not isinstance(value, str) or not value:
        raise ValueError("Each pair must include a non-empty 'aa_seq' or 'sequence'")
    return value


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def write_fasta(path: Path, pairs: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for pair in pairs:
            handle.write(f">{pair_id(pair)}\n{pair_aa(pair)}\n")


def parse_cd_hit_clusters(cluster_path: Path) -> list[list[str]]:
    clusters: list[list[str]] = []
    current: list[str] = []
    with cluster_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = []
                continue
            marker = ">"
            if marker in line:
                current.append(line.split(marker, 1)[1].split("...", 1)[0])
    if current:
        clusters.append(current)
    return clusters


def split_clusters(clusters: list[list[str]], eval_fraction: float) -> tuple[set[str], set[str]]:
    sorted_clusters = sorted(clusters, key=lambda cluster: (-len(cluster), cluster[0] if cluster else ""))
    total = sum(len(cluster) for cluster in sorted_clusters)
    eval_target = max(1, int(total * eval_fraction)) if total > 1 else 0
    train_ids: set[str] = set()
    eval_ids: set[str] = set()

    for cluster in sorted_clusters:
        if len(eval_ids) < eval_target:
            eval_ids.update(cluster)
        else:
            train_ids.update(cluster)

    if not train_ids and eval_ids:
        moved = sorted(eval_ids).pop()
        eval_ids.remove(moved)
        train_ids.add(moved)
    return train_ids, eval_ids


def cd_hit_split(pairs: list[dict[str, Any]], identity: float) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    with tempfile.TemporaryDirectory(prefix="factorforge_run2_") as tmp:
        tmp_dir = Path(tmp)
        input_fasta = tmp_dir / "pairs.fasta"
        output_fasta = tmp_dir / "pairs_cd_hit.fasta"
        write_fasta(input_fasta, pairs)
        subprocess.run(
            [
                "cd-hit",
                "-i",
                str(input_fasta),
                "-o",
                str(output_fasta),
                "-c",
                str(identity),
                "-n",
                "5",
                "-T",
                "0",
                "-M",
                "0",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        clusters = parse_cd_hit_clusters(output_fasta.with_suffix(".fasta.clstr"))

    train_ids, eval_ids = split_clusters(clusters, eval_fraction=0.1)
    by_id = {pair_id(pair): pair for pair in pairs}
    return [by_id[item] for item in sorted(train_ids)], [by_id[item] for item in sorted(eval_ids)]


def fallback_length_split(pairs: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    print("WARNING: cd-hit not found; using length-based fallback split.", file=sys.stderr)
    ordered = sorted(pairs, key=lambda pair: (len(pair_aa(pair)), pair_id(pair)))
    train: list[dict[str, Any]] = []
    eval_rows: list[dict[str, Any]] = []
    for index, pair in enumerate(ordered):
        if index % 10 == 0:
            eval_rows.append(pair)
        else:
            train.append(pair)
    if not train and eval_rows:
        train.append(eval_rows.pop())
    return train, eval_rows


def verify_overlap(train_path: Path, eval_path: Path) -> bool:
    train_ids = {pair_id(pair) for pair in load_pairs(train_path)} if train_path.exists() else set()
    eval_ids = {pair_id(pair) for pair in load_pairs(eval_path)} if eval_path.exists() else set()
    overlap = train_ids & eval_ids
    print(f"train: {len(train_ids)}")
    print(f"eval: {len(eval_ids)}")
    print(f"overlap: {len(overlap)}")
    if overlap:
        print("overlap_ids:", ", ".join(sorted(overlap)[:20]))
    return not overlap


def main() -> int:
    parser = argparse.ArgumentParser(description="Create homology-aware train/eval split.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--train-output", type=Path, default=DEFAULT_TRAIN)
    parser.add_argument("--eval-output", type=Path, default=DEFAULT_EVAL)
    parser.add_argument("--identity", type=float, default=0.7)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    if args.verify and args.train_output.exists() and args.eval_output.exists():
        return 0 if verify_overlap(args.train_output, args.eval_output) else 1

    if not args.input.exists():
        print(f"Input JSONL not found: {args.input}")
        return 0 if args.verify else 1

    pairs = load_pairs(args.input)
    if shutil.which("cd-hit"):
        train_rows, eval_rows = cd_hit_split(pairs, args.identity)
    else:
        train_rows, eval_rows = fallback_length_split(pairs)

    write_jsonl(args.train_output, train_rows)
    write_jsonl(args.eval_output, eval_rows)
    print(f"Saved train: {len(train_rows)} -> {args.train_output}")
    print(f"Saved eval: {len(eval_rows)} -> {args.eval_output}")

    if args.verify:
        return 0 if verify_overlap(args.train_output, args.eval_output) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
