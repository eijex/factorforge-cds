"""Evaluate a Run 4 checkpoint with constrained decoding against v2 pseudo-labels."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "3_training"))

from factorforge.engines.v3.inference.constrained_decoder import constrained_greedy_decode  # noqa: E402
from factorforge.engines.v3.metrics import load_codon_usage_table  # noqa: E402
from factorforge.engines.v3.tokenizer import CodonTokenizer  # noqa: E402
from factorforge.ml.metrics import (  # noqa: E402
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    detect_forbidden_motifs,
    detect_homopolymers,
    detect_repeats,
)
from factorforge.utils.validation import validate_candidate_sequence  # noqa: E402
from train_v3_esm2_bart import build_model, load_config  # noqa: E402


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _protein_sequence(row: dict[str, Any]) -> str:
    sequence = row.get("amino_acid_sequence") or row.get("protein_sequence") or row.get("sequence")
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("Evaluation row is missing amino acid sequence metadata")
    return sequence


def _dna_sequence(row: dict[str, Any]) -> str:
    sequence = row.get("dna_sequence") or row.get("dna_seq")
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("Evaluation row is missing v2 pseudo-label DNA sequence")
    return sequence


def _metric_row(
    protein_id: str,
    engine: str,
    protein: str,
    dna: str,
    cai_weights: dict[str, float],
) -> dict[str, Any]:
    validator = validate_candidate_sequence(protein, dna)
    windows = calculate_gc_windows(dna)
    window_values = [float(window["gc"]) for window in windows]
    return {
        "protein_id": protein_id,
        "engine": engine,
        "cai": calculate_cai(dna, cai_weights),
        "gc_global": calculate_gc(dna),
        "local_gc_min": min(window_values) if window_values else 0.0,
        "local_gc_max": max(window_values) if window_values else 0.0,
        "amino_acid_identity": validator["amino_acid_identity"],
        "validator_pass": validator["passed"],
        "internal_stop_count": validator["internal_stop_count"],
        "invalid_codon_count": validator["invalid_codon_count"],
        "forbidden_motif_count": len(detect_forbidden_motifs(dna, [])),
        "homopolymer_count": len(detect_homopolymers(dna)),
        "repeat_count": len(detect_repeats(dna)),
    }


def evaluate(
    config_path: Path,
    checkpoint_path: Path,
    eval_file: Path | None,
    embeddings_dir: Path | None,
    output_dir: Path,
    limit: int | None = None,
    progress_every: int = 25,
) -> dict[str, Any]:
    cfg = load_config(str(config_path))
    eval_path = eval_file or Path(cfg["data"]["eval_file"])
    embedding_path = embeddings_dir or Path(cfg["paths"]["embeddings_dir"])
    rows = _read_jsonl(eval_path)
    if limit is not None:
        rows = rows[:limit]
    if not rows:
        raise ValueError("No evaluation rows available")

    tokenizer = CodonTokenizer.default()
    table = load_codon_usage_table()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = build_model(cfg).to(device)
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
    model.eval()

    output_dir.mkdir(parents=True, exist_ok=True)
    comparison_rows: list[dict[str, Any]] = []
    v3_metrics: list[dict[str, Any]] = []
    start_time = time.time()
    print(
        f"Run 4 constrained eval: {len(rows)} sequence(s) on {device}; "
        "decoding is autoregressive and may take a while.",
        flush=True,
    )

    for index, row in enumerate(rows, start=1):
        protein_id = str(row["protein_id"])
        protein = _protein_sequence(row)
        embedding_file = embedding_path / f"{protein_id}.pt"
        if not embedding_file.exists():
            raise FileNotFoundError(f"Missing embedding: {embedding_file}")
        embedding = torch.load(embedding_file, map_location=device, weights_only=True)["embeddings"]
        with torch.no_grad():
            token_ids = constrained_greedy_decode(
                model,
                embedding.to(device).unsqueeze(0),
                protein,
                tokenizer,
            )
        v3_dna = tokenizer.decode(token_ids.squeeze(0).tolist(), skip_special_tokens=True)
        v2 = _metric_row(protein_id, "v2", protein, _dna_sequence(row), table.codon_weights)
        v3 = _metric_row(protein_id, "v3", protein, v3_dna, table.codon_weights)
        v3_metrics.append(v3)
        comparison_rows.append(
            {
                "protein_id": protein_id,
                "v2_cai": f"{v2['cai']:.8f}",
                "v3_cai": f"{v3['cai']:.8f}",
                "v3_cai_delta_pct": (
                    f"{((v3['cai'] - v2['cai']) / v2['cai']) * 100:.4f}"
                    if v2["cai"]
                    else "0.0000"
                ),
                "v2_gc_global": f"{v2['gc_global']:.8f}",
                "v3_gc_global": f"{v3['gc_global']:.8f}",
                "v2_local_gc_min": f"{v2['local_gc_min']:.8f}",
                "v2_local_gc_max": f"{v2['local_gc_max']:.8f}",
                "v3_local_gc_min": f"{v3['local_gc_min']:.8f}",
                "v3_local_gc_max": f"{v3['local_gc_max']:.8f}",
                "v3_amino_acid_identity": f"{v3['amino_acid_identity']:.8f}",
                "v3_validator_pass": v3["validator_pass"],
                "v3_internal_stop_count": v3["internal_stop_count"],
                "v3_invalid_codon_count": v3["invalid_codon_count"],
                "v3_forbidden_motif_count": v3["forbidden_motif_count"],
                "v3_homopolymer_count": v3["homopolymer_count"],
                "v3_repeat_count": v3["repeat_count"],
            }
        )
        if progress_every > 0 and (index == 1 or index % progress_every == 0 or index == len(rows)):
            elapsed = time.time() - start_time
            rate = index / elapsed if elapsed > 0 else 0.0
            remaining = (len(rows) - index) / rate if rate > 0 else 0.0
            print(
                f"Decoded {index}/{len(rows)} | elapsed {elapsed / 60:.1f} min | "
                f"eta {remaining / 60:.1f} min",
                flush=True,
            )

    comparison_path = output_dir / "alpha_run1_comparison.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0].keys()))
        writer.writeheader()
        writer.writerows(comparison_rows)

    summary = {
        "decoded_count": len(v3_metrics),
        "validator_pass_rate": sum(1 for item in v3_metrics if item["validator_pass"])
        / len(v3_metrics),
        "min_amino_acid_identity": min(float(item["amino_acid_identity"]) for item in v3_metrics),
        "min_gc_global": min(float(item["gc_global"]) for item in v3_metrics),
        "max_gc_global": max(float(item["gc_global"]) for item in v3_metrics),
        "min_v3_to_v2_cai_ratio": min(
            (
                float(row["v3_cai"]) / float(row["v2_cai"])
                if float(row["v2_cai"]) > 0
                else 1.0
            )
            for row in comparison_rows
        ),
        "comparison_csv": str(comparison_path),
    }
    (output_dir / "alpha_run1_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run constrained Run 4 evaluation.")
    parser.add_argument("--config", type=Path, default=ROOT / "configs" / "v3_training_config_alpha_run1.yml")
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--eval-file", type=Path, default=None)
    parser.add_argument("--embeddings-dir", type=Path, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "experiments" / "results" / "alpha_run1",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=25)
    args = parser.parse_args()
    print(
        json.dumps(
            evaluate(
                config_path=args.config,
                checkpoint_path=args.checkpoint,
                eval_file=args.eval_file,
                embeddings_dir=args.embeddings_dir,
                output_dir=args.output_dir,
                limit=args.limit,
                progress_every=args.progress_every,
            ),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
