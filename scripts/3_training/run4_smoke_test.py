"""Run the FactorForge v3-alpha Run 4 smoke test.

This runner is intentionally smoke-only: it creates a small synthetic protein
panel, generates v2 pseudo-labels, trains a small constrained v3 decoder for a
bounded number of steps, then writes machine-readable v2/v3 comparison metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from pathlib import Path
from typing import Any

import torch
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "1_data_preparation"))
sys.path.insert(0, str(ROOT / "scripts" / "3_training"))

from generate_v2_pseudolabels import generate_v2_pseudolabels  # noqa: E402
from factorforge.engines.v3.inference.constrained_decoder import (  # noqa: E402
    constrained_greedy_decode,
)
from factorforge.engines.v3.metrics import load_codon_usage_table  # noqa: E402
from factorforge.engines.v3.modeling_bart_decoder import BartDecoderSkeleton  # noqa: E402
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
from dataset import CodonDataset, collate_fn  # noqa: E402
from train_v3_esm2_bart import (  # noqa: E402
    _codon_gc_vector,
    _codon_log_cai_vector,
    _codon_token_mask,
    _compute_bounded_gc_penalty,
    _compute_ce_loss,
    _compute_expected_gc,
    _compute_expected_log_cai,
)


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
GC_STABLE_AMINO_ACIDS = "AGPW"


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _write_fasta(path: Path, count: int, seed: int) -> None:
    if not 50 <= count <= 200:
        raise ValueError("Run 4 smoke requires 50-200 protein sequences")

    rng = random.Random(seed)
    records: list[tuple[str, str]] = []
    while len(records) < count:
        length = rng.randint(12, 48)
        sequence = "M" + "".join(
            rng.choice(GC_STABLE_AMINO_ACIDS) for _ in range(length - 1)
        )
        records.append((f"gc_stable_synthetic_{len(records) + 1:03d}", sequence))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for protein_id, sequence in records[:count]:
            handle.write(f">{protein_id}\n{sequence}\n")


def _embedding_for(protein_sequence: str, dim: int) -> torch.Tensor:
    tensor = torch.zeros((len(protein_sequence), dim), dtype=torch.float32)
    for pos, aa in enumerate(protein_sequence):
        aa_index = AMINO_ACIDS.index(aa) if aa in AMINO_ACIDS else 0
        tensor[pos, aa_index % dim] = 1.0
        if dim > 20:
            tensor[pos, 20] = pos / max(len(protein_sequence) - 1, 1)
        if dim > 21:
            tensor[pos, 21] = aa_index / (len(AMINO_ACIDS) - 1)
    return tensor


def _write_embeddings(rows: list[dict[str, Any]], embeddings_dir: Path, dim: int) -> None:
    embeddings_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        protein = str(row["protein_sequence"])
        torch.save(
            {"embeddings": _embedding_for(protein, dim)},
            embeddings_dir / f"{row['protein_id']}.pt",
        )


def _smoke_cfg(cfg: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    smoke = json.loads(json.dumps(cfg))
    smoke["paths"]["embeddings_dir"] = str(args.work_dir / "embeddings")
    smoke["data"]["train_file"] = str(args.work_dir / "run4_v2_pseudolabels_train.jsonl")
    smoke["data"]["eval_file"] = str(args.work_dir / "run4_v2_pseudolabels_eval.jsonl")
    smoke["training"]["max_steps"] = args.max_steps
    smoke["training"]["batch_size"] = args.batch_size
    smoke["training"]["eval_every"] = args.log_every
    smoke["training"]["checkpoint_every"] = args.max_steps + 1
    smoke["training"]["early_stopping_patience"] = 0

    # CPU smoke override: keeps Run 4 loss/config semantics but avoids full-training scale.
    smoke["bart"] = {
        **smoke["bart"],
        "d_model": args.d_model,
        "decoder_layers": args.decoder_layers,
        "decoder_attention_heads": args.decoder_attention_heads,
        "ffn_dim": args.ffn_dim,
        "max_position_embeddings": 128,
        "dropout": 0.10,
    }
    return smoke


def _build_model(cfg: dict[str, Any]) -> BartDecoderSkeleton:
    tokenizer = CodonTokenizer.default()
    bart_cfg = cfg["bart"]
    return BartDecoderSkeleton(
        vocab_size=len(tokenizer.token_to_id),
        d_model=int(bart_cfg["d_model"]),
        encoder_dim=int(bart_cfg["encoder_dim"]),
        decoder_layers=int(bart_cfg["decoder_layers"]),
        decoder_attention_heads=int(bart_cfg["decoder_attention_heads"]),
        ffn_dim=int(bart_cfg["ffn_dim"]),
        max_position_embeddings=int(bart_cfg["max_position_embeddings"]),
        dropout=float(bart_cfg["dropout"]),
        pad_token_id=tokenizer.pad_token_id,
        bos_token_id=tokenizer.bos_token_id,
        eos_token_id=tokenizer.eos_token_id,
    )


def _synonym_mask_coverage(batch: dict[str, torch.Tensor]) -> float:
    labels = batch["labels"]
    valid_positions = labels.ne(-100) & labels.ge(5)
    if not bool(valid_positions.any()):
        return 0.0
    valid_mask = batch["synonym_mask"].any(dim=-1)
    return float((valid_mask & valid_positions).sum().item() / valid_positions.sum().item())


def _train_smoke(
    cfg: dict[str, Any],
    loss_log_path: Path,
    *,
    seed: int,
) -> tuple[BartDecoderSkeleton, int]:
    torch.manual_seed(seed)
    device = torch.device("cpu")
    tokenizer = CodonTokenizer.default()
    model = _build_model(cfg).to(device)
    codon_gc = _codon_gc_vector(tokenizer, device)
    codon_mask = _codon_token_mask(tokenizer, device)
    codon_log_cai = _codon_log_cai_vector(tokenizer, device)
    loss_cfg = cfg["loss"]
    training_cfg = cfg["training"]

    dataset = CodonDataset(
        training_jsonl=cfg["data"]["train_file"],
        embeddings_dir=cfg["paths"]["embeddings_dir"],
        codon_to_id=tokenizer.token_to_id,
        bos_token_id=tokenizer.bos_token_id,
        eos_token_id=tokenizer.eos_token_id,
        unk_token_id=tokenizer.unk_token_id,
        max_length=int(cfg["bart"]["max_position_embeddings"]) - 2,
    )
    generator = torch.Generator().manual_seed(seed)
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=int(training_cfg["batch_size"]),
        shuffle=True,
        collate_fn=collate_fn,
        generator=generator,
        num_workers=0,
    )
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(training_cfg["learning_rate"]),
        weight_decay=float(training_cfg.get("weight_decay", 0.0)),
    )

    max_steps = int(training_cfg["max_steps"])
    loss_log_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "step",
        "CE_to_pseudo_label",
        "expected_GC",
        "bounded_GC_penalty",
        "expected_log_CAI",
        "total_loss",
        "synonym_mask_coverage",
    ]
    with loss_log_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        step = 0
        while step < max_steps:
            for batch in loader:
                if step >= max_steps:
                    break
                batch = {key: value.to(device) for key, value in batch.items()}
                logits = model(
                    encoder_hidden_states=batch["encoder_hidden_states"],
                    decoder_input_ids=batch["decoder_input_ids"],
                )
                ce = _compute_ce_loss(
                    logits,
                    batch["labels"],
                    label_smoothing=float(training_cfg.get("label_smoothing", 0.0)),
                )
                expected_gc = _compute_expected_gc(
                    logits,
                    batch["labels"],
                    codon_gc,
                    codon_mask,
                    batch["synonym_mask"],
                )
                bounded_gc = _compute_bounded_gc_penalty(
                    expected_gc,
                    gc_low=float(loss_cfg["gc_low"]),
                    gc_high=float(loss_cfg["gc_high"]),
                    lambda_low=float(loss_cfg.get("gc_lambda_low", 1.0)),
                    lambda_high=float(loss_cfg.get("gc_lambda_high", 1.0)),
                )
                expected_log_cai = _compute_expected_log_cai(
                    logits,
                    batch["labels"],
                    codon_log_cai,
                    codon_mask,
                    batch["synonym_mask"],
                )
                total = (
                    float(loss_cfg.get("ce_weight", 1.0)) * ce
                    + float(loss_cfg.get("gc_weight", 0.0)) * bounded_gc
                    - float(loss_cfg.get("expected_log_cai_weight", 0.0)) * expected_log_cai
                )

                optimizer.zero_grad()
                total.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

                step += 1
                writer.writerow(
                    {
                        "step": step,
                        "CE_to_pseudo_label": f"{ce.item():.8f}",
                        "expected_GC": f"{expected_gc.item():.8f}",
                        "bounded_GC_penalty": f"{bounded_gc.item():.8f}",
                        "expected_log_CAI": f"{expected_log_cai.item():.8f}",
                        "total_loss": f"{total.item():.8f}",
                        "synonym_mask_coverage": f"{_synonym_mask_coverage(batch):.8f}",
                    }
                )
    return model, max_steps


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


def _decode_and_compare(
    model: BartDecoderSkeleton,
    rows: list[dict[str, Any]],
    embeddings_dir: Path,
    comparison_path: Path,
) -> dict[str, Any]:
    tokenizer = CodonTokenizer.default()
    table = load_codon_usage_table()
    model.eval()
    paired_rows: list[dict[str, Any]] = []
    v3_rows: list[dict[str, Any]] = []

    for row in rows:
        protein_id = str(row["protein_id"])
        protein = str(row["protein_sequence"])
        emb = torch.load(embeddings_dir / f"{protein_id}.pt", map_location="cpu", weights_only=True)
        with torch.no_grad():
            token_ids = constrained_greedy_decode(
                model,
                emb["embeddings"].unsqueeze(0),
                protein,
                tokenizer,
            )
        v3_dna = tokenizer.decode(token_ids.squeeze(0).tolist(), skip_special_tokens=True)
        v2_metrics = _metric_row(protein_id, "v2", protein, str(row["dna_seq"]), table.codon_weights)
        v3_metrics = _metric_row(protein_id, "v3", protein, v3_dna, table.codon_weights)
        v3_rows.append(v3_metrics)
        paired_rows.append(
            {
                "protein_id": protein_id,
                "v2_cai": f"{v2_metrics['cai']:.8f}",
                "v3_cai": f"{v3_metrics['cai']:.8f}",
                "v3_cai_delta_pct": (
                    f"{((v3_metrics['cai'] - v2_metrics['cai']) / v2_metrics['cai']) * 100:.4f}"
                    if v2_metrics["cai"]
                    else "0.0000"
                ),
                "v2_gc_global": f"{v2_metrics['gc_global']:.8f}",
                "v3_gc_global": f"{v3_metrics['gc_global']:.8f}",
                "v2_local_gc_min": f"{v2_metrics['local_gc_min']:.8f}",
                "v2_local_gc_max": f"{v2_metrics['local_gc_max']:.8f}",
                "v3_local_gc_min": f"{v3_metrics['local_gc_min']:.8f}",
                "v3_local_gc_max": f"{v3_metrics['local_gc_max']:.8f}",
                "v3_amino_acid_identity": f"{v3_metrics['amino_acid_identity']:.8f}",
                "v3_validator_pass": v3_metrics["validator_pass"],
                "v3_internal_stop_count": v3_metrics["internal_stop_count"],
                "v3_invalid_codon_count": v3_metrics["invalid_codon_count"],
                "v3_forbidden_motif_count": v3_metrics["forbidden_motif_count"],
                "v3_homopolymer_count": v3_metrics["homopolymer_count"],
                "v3_repeat_count": v3_metrics["repeat_count"],
            }
        )

    comparison_path.parent.mkdir(parents=True, exist_ok=True)
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(paired_rows[0].keys()))
        writer.writeheader()
        writer.writerows(paired_rows)

    validator_pass_rate = sum(1 for row in v3_rows if row["validator_pass"]) / len(v3_rows)
    min_gc = min(float(row["gc_global"]) for row in v3_rows)
    min_identity = min(float(row["amino_acid_identity"]) for row in v3_rows)
    min_cai_ratio = min(
        (
            float(paired["v3_cai"]) / float(paired["v2_cai"])
            if float(paired["v2_cai"]) > 0
            else 1.0
        )
        for paired in paired_rows
    )
    return {
        "decoded_count": len(v3_rows),
        "validator_pass_rate": validator_pass_rate,
        "min_gc_global": min_gc,
        "min_amino_acid_identity": min_identity,
        "min_v3_to_v2_cai_ratio": min_cai_ratio,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    cfg = _smoke_cfg(_load_yaml(args.config), args)

    fasta_path = args.work_dir / "run4_smoke_proteins.fasta"
    all_jsonl = args.work_dir / "run4_v2_pseudolabels.jsonl"
    _write_fasta(fasta_path, args.sequence_count, args.seed)
    counts = generate_v2_pseudolabels(
        fasta_path=fasta_path,
        output_path=all_jsonl,
        profile=str(cfg["data"]["pseudo_label_profile"]),
        max_length=int(cfg["data"]["max_protein_length"]),
        require_validator_pass=True,
    )
    rows = _read_jsonl(all_jsonl)
    if len(rows) < 50:
        raise RuntimeError(f"Expected at least 50 pseudo-label rows, got {len(rows)}")

    split = max(1, int(len(rows) * 0.8))
    train_rows = rows[:split]
    eval_rows = rows[split:] or rows[-1:]
    _write_jsonl(Path(cfg["data"]["train_file"]), train_rows)
    _write_jsonl(Path(cfg["data"]["eval_file"]), eval_rows)
    _write_embeddings(rows, Path(cfg["paths"]["embeddings_dir"]), int(cfg["bart"]["encoder_dim"]))

    model, steps = _train_smoke(cfg, args.results_dir / "run4_smoke_loss_log.csv", seed=args.seed)
    comparison = _decode_and_compare(
        model,
        eval_rows,
        Path(cfg["paths"]["embeddings_dir"]),
        args.results_dir / "run4_smoke_comparison.csv",
    )
    criteria = {
        "smoke_training_completed": True,
        "validator_pass_rate_100pct": comparison["validator_pass_rate"] == 1.0,
        "no_gc_collapse_gc_ge_40": comparison["min_gc_global"] >= 40.0,
        "amino_acid_identity_100pct": comparison["min_amino_acid_identity"] == 1.0,
        "v3_cai_within_minus_10pct_of_v2": comparison["min_v3_to_v2_cai_ratio"] >= 0.90,
        "metrics_saved": True,
    }
    summary = {
        "job": "025-medium-run4-smoke-test",
        "source_config": str(args.config),
        "results_dir": str(args.results_dir),
        "work_dir": str(args.work_dir),
        "sequence_count_requested": args.sequence_count,
        "pseudo_label_counts": counts,
        "train_rows": len(train_rows),
        "eval_rows": len(eval_rows),
        "training_steps": steps,
        "smoke_model_overrides": {
            "d_model": args.d_model,
            "decoder_layers": args.decoder_layers,
            "decoder_attention_heads": args.decoder_attention_heads,
            "ffn_dim": args.ffn_dim,
            "reason": "CPU smoke run; full Run 4 training remains out of scope",
        },
        "comparison": comparison,
        "acceptance_criteria": criteria,
        "verdict": "pass" if all(criteria.values()) else "fail",
        "full_run4_recommendation": (
            "proceed_to_approval_review" if all(criteria.values()) else "do_not_proceed"
        ),
    }
    (args.results_dir / "run4_smoke_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run bounded Run 4 smoke training and validation.")
    parser.add_argument("--config", type=Path, default=ROOT / "configs" / "v3_training_config_alpha_run1.yml")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=ROOT / "experiments" / "results" / "run4_smoke",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=ROOT / "data" / "training" / "run4_smoke_work",
    )
    parser.add_argument("--sequence-count", type=int, default=60)
    parser.add_argument("--max-steps", type=int, default=500)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--log-every", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--d-model", type=int, default=64)
    parser.add_argument("--decoder-layers", type=int, default=2)
    parser.add_argument("--decoder-attention-heads", type=int, default=4)
    parser.add_argument("--ffn-dim", type=int, default=128)
    args = parser.parse_args()

    if args.max_steps < 500 or args.max_steps > 2000:
        raise ValueError("--max-steps must be between 500 and 2000 for this smoke test")

    summary = run(args)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["verdict"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
