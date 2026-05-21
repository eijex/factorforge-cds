"""Train FactorForge v3 with per-token ESM2 embeddings and a BART decoder."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import torch
import torch.nn.functional as F
import yaml
from torch.optim import AdamW
from torch.utils.data import DataLoader, random_split

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from factorforge.engines.v3.modeling_bart_decoder import BartDecoderSkeleton
from factorforge.engines.v3.metrics import load_codon_usage_table
from factorforge.engines.v3.tokenizer import CodonTokenizer
from dataset import CodonDataset, collate_fn


def load_config(config_path: str) -> dict:
    with open(config_path, encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def build_model(cfg: dict) -> BartDecoderSkeleton:
    codon_tok = CodonTokenizer.default()
    bart_cfg = cfg.get("bart") or cfg["model"]
    return BartDecoderSkeleton(
        vocab_size=len(codon_tok.token_to_id),
        d_model=bart_cfg["d_model"],
        encoder_dim=bart_cfg["encoder_dim"],
        decoder_layers=bart_cfg["decoder_layers"],
        decoder_attention_heads=bart_cfg["decoder_attention_heads"],
        ffn_dim=bart_cfg["ffn_dim"],
        max_position_embeddings=bart_cfg["max_position_embeddings"],
        dropout=bart_cfg["dropout"],
        pad_token_id=codon_tok.pad_token_id,
        bos_token_id=codon_tok.bos_token_id,
        eos_token_id=codon_tok.eos_token_id,
    )


def _compute_ce_loss(
    logits: torch.Tensor,
    labels: torch.Tensor,
    label_smoothing: float = 0.0,
) -> torch.Tensor:
    return F.cross_entropy(
        logits.reshape(-1, logits.size(-1)),
        labels.reshape(-1),
        ignore_index=-100,
        label_smoothing=label_smoothing,
    )


def _codon_gc_vector(codon_tok: CodonTokenizer, device: torch.device) -> torch.Tensor:
    values: list[float] = []
    for idx in range(len(codon_tok.id_to_token)):
        token = codon_tok.id_to_token[idx]
        if token in codon_tok.special_tokens or len(token) != 3:
            values.append(0.0)
        else:
            values.append((token.count("G") + token.count("C")) / 3.0)
    return torch.tensor(values, dtype=torch.float32, device=device)


def _codon_token_mask(codon_tok: CodonTokenizer, device: torch.device) -> torch.Tensor:
    values: list[bool] = []
    for idx in range(len(codon_tok.id_to_token)):
        token = codon_tok.id_to_token[idx]
        values.append(token not in codon_tok.special_tokens and len(token) == 3)
    return torch.tensor(values, dtype=torch.bool, device=device)


def _codon_log_cai_vector(
    codon_tok: CodonTokenizer,
    device: torch.device,
    epsilon: float = 1e-8,
) -> torch.Tensor:
    table = load_codon_usage_table()
    values: list[float] = []
    for idx in range(len(codon_tok.id_to_token)):
        token = codon_tok.id_to_token[idx]
        if token in codon_tok.special_tokens or len(token) != 3:
            values.append(0.0)
        else:
            values.append(float(torch.log(torch.tensor(max(table.codon_weights.get(token, 0.0), epsilon)))))
    return torch.tensor(values, dtype=torch.float32, device=device)


def _valid_label_positions(labels: torch.Tensor) -> torch.Tensor:
    # Special token IDs are fixed before codon IDs in CodonTokenizer.
    return labels.ne(-100) & labels.ge(5)


def _masked_codon_probs(
    logits: torch.Tensor,
    codon_token_mask: torch.Tensor,
    synonym_mask: torch.Tensor | None = None,
) -> torch.Tensor:
    allowed = codon_token_mask.view(1, 1, -1).expand_as(logits)
    if synonym_mask is not None:
        allowed = allowed & synonym_mask.to(dtype=torch.bool, device=logits.device)
        fallback_allowed = codon_token_mask.view(1, 1, -1).expand_as(logits)
        allowed = torch.where(allowed.any(dim=-1, keepdim=True), allowed, fallback_allowed)
    masked_logits = logits.masked_fill(~allowed, -1.0e9)
    return torch.softmax(masked_logits, dim=-1)


def _compute_expected_gc(
    logits: torch.Tensor,
    labels: torch.Tensor,
    codon_gc: torch.Tensor,
    codon_token_mask: torch.Tensor | None = None,
    synonym_mask: torch.Tensor | None = None,
) -> torch.Tensor:
    valid_positions = _valid_label_positions(labels)
    if not bool(valid_positions.any()):
        return logits.new_tensor(0.0)

    if codon_token_mask is None:
        codon_token_mask = torch.ones(logits.shape[-1], dtype=torch.bool, device=logits.device)
    probs = _masked_codon_probs(logits, codon_token_mask, synonym_mask=synonym_mask)
    expected_gc_by_position = torch.matmul(probs, codon_gc)
    return expected_gc_by_position[valid_positions].mean()


def _compute_expected_log_cai(
    logits: torch.Tensor,
    labels: torch.Tensor,
    codon_log_cai: torch.Tensor,
    codon_token_mask: torch.Tensor | None = None,
    synonym_mask: torch.Tensor | None = None,
) -> torch.Tensor:
    valid_positions = _valid_label_positions(labels)
    if not bool(valid_positions.any()):
        return logits.new_tensor(0.0)

    if codon_token_mask is None:
        codon_token_mask = torch.ones(logits.shape[-1], dtype=torch.bool, device=logits.device)
    probs = _masked_codon_probs(logits, codon_token_mask, synonym_mask=synonym_mask)
    expected_log_cai_by_position = torch.matmul(probs, codon_log_cai)
    return expected_log_cai_by_position[valid_positions].mean()


def _compute_bounded_gc_penalty(
    expected_gc: torch.Tensor,
    gc_low: float = 0.40,
    gc_high: float = 0.55,
    lambda_low: float = 1.0,
    lambda_high: float = 1.0,
) -> torch.Tensor:
    return lambda_low * torch.relu(expected_gc.new_tensor(gc_low) - expected_gc) + (
        lambda_high * torch.relu(expected_gc - expected_gc.new_tensor(gc_high))
    )


def _compute_gc_penalty(
    logits: torch.Tensor,
    labels: torch.Tensor,
    codon_gc: torch.Tensor,
    gc_target: float,
    codon_token_mask: torch.Tensor | None = None,
) -> torch.Tensor:
    """Backward-compatible single-target GC penalty."""
    expected_gc = _compute_expected_gc(logits, labels, codon_gc, codon_token_mask)
    return torch.abs(expected_gc - gc_target)


def _compute_loss(
    logits: torch.Tensor,
    labels: torch.Tensor,
    codon_gc: torch.Tensor,
    loss_cfg: dict,
    training_cfg: dict,
    synonym_mask: torch.Tensor | None = None,
) -> torch.Tensor:
    ce_weight = float(loss_cfg.get("ce_weight", 1.0))
    gc_weight = float(loss_cfg.get("gc_weight", 0.0))
    gc_low = float(loss_cfg.get("gc_low", 0.40))
    gc_high = float(loss_cfg.get("gc_high", 0.55))
    lambda_low = float(loss_cfg.get("gc_lambda_low", 1.0))
    lambda_high = float(loss_cfg.get("gc_lambda_high", 1.0))
    cai_weight = float(loss_cfg.get("cai_weight", loss_cfg.get("expected_log_cai_weight", 0.0)))
    label_smoothing = float(training_cfg.get("label_smoothing", 0.0))

    ce_loss = _compute_ce_loss(logits, labels, label_smoothing=label_smoothing)
    loss = ce_weight * ce_loss
    codon_tok = CodonTokenizer.default()
    codon_mask = _codon_token_mask(codon_tok, logits.device)
    if gc_weight > 0:
        expected_gc = _compute_expected_gc(logits, labels, codon_gc, codon_mask, synonym_mask)
        gc_loss = _compute_bounded_gc_penalty(
            expected_gc,
            gc_low=gc_low,
            gc_high=gc_high,
            lambda_low=lambda_low,
            lambda_high=lambda_high,
        )
        loss = loss + gc_weight * gc_loss
    if cai_weight > 0:
        codon_log_cai = _codon_log_cai_vector(codon_tok, logits.device)
        expected_log_cai = _compute_expected_log_cai(
            logits,
            labels,
            codon_log_cai,
            codon_mask,
            synonym_mask,
        )
        loss = loss - cai_weight * expected_log_cai
    return loss


def _compute_loss_components(
    logits: torch.Tensor,
    labels: torch.Tensor,
    codon_gc: torch.Tensor,
    loss_cfg: dict,
    training_cfg: dict,
    synonym_mask: torch.Tensor | None = None,
) -> dict[str, torch.Tensor]:
    codon_tok = CodonTokenizer.default()
    codon_mask = _codon_token_mask(codon_tok, logits.device)
    ce = _compute_ce_loss(
        logits,
        labels,
        label_smoothing=float(training_cfg.get("label_smoothing", 0.0)),
    )
    expected_gc = _compute_expected_gc(logits, labels, codon_gc, codon_mask, synonym_mask)
    bounded_gc = _compute_bounded_gc_penalty(
        expected_gc,
        gc_low=float(loss_cfg.get("gc_low", 0.40)),
        gc_high=float(loss_cfg.get("gc_high", 0.55)),
        lambda_low=float(loss_cfg.get("gc_lambda_low", 1.0)),
        lambda_high=float(loss_cfg.get("gc_lambda_high", 1.0)),
    )
    expected_log_cai = _compute_expected_log_cai(
        logits,
        labels,
        _codon_log_cai_vector(codon_tok, logits.device),
        codon_mask,
        synonym_mask,
    )
    total = (
        float(loss_cfg.get("ce_weight", 1.0)) * ce
        + float(loss_cfg.get("gc_weight", 0.0)) * bounded_gc
        - float(loss_cfg.get("cai_weight", loss_cfg.get("expected_log_cai_weight", 0.0)))
        * expected_log_cai
    )
    return {
        "CE_to_pseudo_label": ce,
        "expected_GC": expected_gc,
        "bounded_GC_penalty": bounded_gc,
        "expected_log_CAI": expected_log_cai,
        "total_loss": total,
    }


def _synonym_mask_coverage(labels: torch.Tensor, synonym_mask: torch.Tensor | None) -> float:
    valid_positions = _valid_label_positions(labels)
    if synonym_mask is None or not bool(valid_positions.any()):
        return 0.0
    valid_mask = synonym_mask.any(dim=-1)
    return float((valid_mask & valid_positions).sum().item() / valid_positions.sum().item())


def _training_file(cfg: dict) -> str:
    return cfg.get("data", {}).get("train_file") or cfg["paths"]["training_data"]


def _eval_file(cfg: dict) -> str | None:
    return cfg.get("data", {}).get("eval_file")


def _bart_cfg(cfg: dict) -> dict:
    return cfg.get("bart") or cfg["model"]


def train(cfg: dict, dry_run: bool = False, loss_log: str | None = None) -> None:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    checkpoint_dir = Path(cfg["paths"]["checkpoint_dir"])
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    model = build_model(cfg).to(device)
    total_params = sum(param.numel() for param in model.parameters())
    print(f"Model parameters: {total_params:,}")

    if dry_run:
        print("Dry-run: model initialized successfully. Exiting.")
        return

    codon_tok = CodonTokenizer.default()
    codon_gc = _codon_gc_vector(codon_tok, device)
    bart_cfg = _bart_cfg(cfg)
    training_cfg = cfg["training"]
    loss_cfg = cfg.get("loss", {})
    dataset = CodonDataset(
        training_jsonl=_training_file(cfg),
        embeddings_dir=cfg["paths"]["embeddings_dir"],
        codon_to_id=codon_tok.token_to_id,
        bos_token_id=codon_tok.bos_token_id,
        eos_token_id=codon_tok.eos_token_id,
        unk_token_id=codon_tok.unk_token_id,
        max_length=bart_cfg["max_position_embeddings"] - 2,
    )
    if len(dataset) < 2:
        raise ValueError("Training requires at least 2 samples after embedding matching")

    eval_file = _eval_file(cfg)
    if eval_file and Path(eval_file).exists():
        train_ds = dataset
        val_ds = CodonDataset(
            training_jsonl=eval_file,
            embeddings_dir=cfg["paths"]["embeddings_dir"],
            codon_to_id=codon_tok.token_to_id,
            bos_token_id=codon_tok.bos_token_id,
            eos_token_id=codon_tok.eos_token_id,
            unk_token_id=codon_tok.unk_token_id,
            max_length=bart_cfg["max_position_embeddings"] - 2,
        )
    else:
        train_size = int(len(dataset) * cfg["data"]["train_split"])
        train_size = min(max(train_size, 1), len(dataset) - 1)
        val_size = len(dataset) - train_size
        generator = torch.Generator().manual_seed(training_cfg["seed"])
        train_ds, val_ds = random_split(dataset, [train_size, val_size], generator=generator)

    train_loader = DataLoader(
        train_ds,
        batch_size=training_cfg["batch_size"],
        shuffle=True,
        collate_fn=collate_fn,
        num_workers=2,
    )
    val_loader = DataLoader(
        val_ds,
        batch_size=training_cfg["batch_size"],
        collate_fn=collate_fn,
        num_workers=2,
    )

    optimizer = AdamW(
        model.parameters(),
        lr=training_cfg["learning_rate"],
        weight_decay=float(training_cfg.get("weight_decay", 0.0)),
    )

    start_step = 0
    latest_ckpt = sorted(checkpoint_dir.glob("step_*.pt"))
    if latest_ckpt:
        checkpoint = torch.load(latest_ckpt[-1], map_location=device, weights_only=True)
        model.load_state_dict(checkpoint["model"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        start_step = checkpoint["step"]
        print(f"Resumed from step {start_step}")

    model.train()
    step = start_step
    max_steps = training_cfg["max_steps"]
    best_eval_loss = float("inf")
    stale_evals = 0
    early_stopping_patience = training_cfg.get("early_stopping_patience")
    loss_log_handle = None
    loss_writer = None
    if loss_log:
        loss_log_path = Path(loss_log)
        loss_log_path.parent.mkdir(parents=True, exist_ok=True)
        loss_log_handle = loss_log_path.open("w", newline="", encoding="utf-8")
        loss_writer = csv.DictWriter(
            loss_log_handle,
            fieldnames=[
                "step",
                "split",
                "CE_to_pseudo_label",
                "expected_GC",
                "bounded_GC_penalty",
                "expected_log_CAI",
                "total_loss",
                "synonym_mask_coverage",
            ],
        )
        loss_writer.writeheader()

    try:
        for _epoch in range(1000):
            for batch in train_loader:
                if step >= max_steps:
                    break

                encoder_hidden = batch["encoder_hidden_states"].to(device)
                decoder_ids = batch["decoder_input_ids"].to(device)
                labels = batch["labels"].to(device)
                synonym_mask = batch.get("synonym_mask")
                synonym_mask = synonym_mask.to(device) if synonym_mask is not None else None

                logits = model(
                    encoder_hidden_states=encoder_hidden,
                    decoder_input_ids=decoder_ids,
                )
                components = _compute_loss_components(
                    logits,
                    labels,
                    codon_gc,
                    loss_cfg,
                    training_cfg,
                    synonym_mask=synonym_mask,
                )
                loss = components["total_loss"]

                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

                step += 1

                if step % 100 == 0:
                    print(f"Step {step}/{max_steps} | Loss: {loss.item():.4f}")

                if step % training_cfg["eval_every"] == 0:
                    eval_loss = _evaluate(
                        model,
                        val_loader,
                        device,
                        codon_gc,
                        loss_cfg,
                        training_cfg,
                        loss_writer=loss_writer,
                        step=step,
                    )
                    if loss_log_handle is not None:
                        loss_log_handle.flush()
                    if eval_loss is not None and early_stopping_patience:
                        if eval_loss < best_eval_loss:
                            best_eval_loss = eval_loss
                            stale_evals = 0
                        else:
                            stale_evals += 1
                            if stale_evals >= int(early_stopping_patience):
                                print(f"Early stopping at step {step}; best eval loss: {best_eval_loss:.4f}")
                                step = max_steps
                                break

                if step % training_cfg.get("checkpoint_every", 5000) == 0:
                    ckpt_path = checkpoint_dir / f"step_{step:07d}.pt"
                    torch.save(
                        {
                            "step": step,
                            "model": model.state_dict(),
                            "optimizer": optimizer.state_dict(),
                        },
                        ckpt_path,
                    )
                    print(f"Checkpoint saved: {ckpt_path}")

            if step >= max_steps:
                break
    finally:
        if loss_log_handle is not None:
            loss_log_handle.close()

    final_path = Path(cfg["paths"]["model_output"])
    final_path.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), final_path / "pytorch_model.pt")
    print(f"Training complete. Model saved to {final_path}")


def _evaluate(
    model: BartDecoderSkeleton,
    val_loader: DataLoader,
    device: torch.device,
    codon_gc: torch.Tensor,
    loss_cfg: dict,
    training_cfg: dict,
    loss_writer: csv.DictWriter | None = None,
    step: int | None = None,
) -> float | None:
    model.eval()
    losses: list[float] = []
    with torch.no_grad():
        for batch in val_loader:
            logits = model(
                encoder_hidden_states=batch["encoder_hidden_states"].to(device),
                decoder_input_ids=batch["decoder_input_ids"].to(device),
            )
            synonym_mask = batch["synonym_mask"].to(device) if "synonym_mask" in batch else None
            labels = batch["labels"].to(device)
            components = _compute_loss_components(
                logits,
                labels,
                codon_gc,
                loss_cfg,
                training_cfg,
                synonym_mask=synonym_mask,
            )
            loss = components["total_loss"]
            losses.append(float(loss.item()))
            if loss_writer is not None and step is not None:
                loss_writer.writerow(
                    {
                        "step": step,
                        "split": "eval",
                        "CE_to_pseudo_label": f"{components['CE_to_pseudo_label'].item():.8f}",
                        "expected_GC": f"{components['expected_GC'].item():.8f}",
                        "bounded_GC_penalty": f"{components['bounded_GC_penalty'].item():.8f}",
                        "expected_log_CAI": f"{components['expected_log_CAI'].item():.8f}",
                        "total_loss": f"{components['total_loss'].item():.8f}",
                        "synonym_mask_coverage": f"{_synonym_mask_coverage(labels, synonym_mask):.8f}",
                    }
                )
    if losses:
        mean_loss = sum(losses) / len(losses)
        print(f"Eval loss: {mean_loss:.4f}")
        model.train()
        return mean_loss
    model.train()
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--loss-log", default=None, help="Optional CSV path for eval loss components.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    train(cfg, dry_run=args.dry_run, loss_log=args.loss_log)


if __name__ == "__main__":
    main()
