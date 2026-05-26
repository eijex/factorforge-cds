#!/usr/bin/env python3
"""
Fine-tune a pretrained BART model with full-sequence masking.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import torch
import yaml
from Bio import SeqIO
from torch.utils.data import Dataset
from transformers import (
    BartConfig,
    BartForConditionalGeneration,
    EarlyStoppingCallback,
    Trainer,
    TrainingArguments,
)
from transformers.trainer_callback import TrainerCallback

SPECIAL_TOKENS = ["[PAD]", "[UNK]", "[MASK]", "[START]", "[END]"]


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


class CodonTokenizer:
    def __init__(self, token_map: Dict[str, int]):
        self.token_to_id = token_map
        self.id_to_token = {idx: token for token, idx in token_map.items()}
        self.pad_token_id = token_map["[PAD]"]
        self.unk_token_id = token_map["[UNK]"]
        self.mask_token_id = token_map["[MASK]"]
        self.start_token_id = token_map["[START]"]
        self.end_token_id = token_map["[END]"]
        self.codon_ids = [
            idx
            for token, idx in token_map.items()
            if len(token) == 3 and all(base in "ACGT" for base in token)
        ]

    @classmethod
    def from_json(cls, path: Path) -> "CodonTokenizer":
        with path.open("r", encoding="utf-8") as handle:
            token_map = json.load(handle)
        missing = [token for token in SPECIAL_TOKENS if token not in token_map]
        if missing:
            raise ValueError(f"Tokenizer missing special tokens: {missing}")
        return cls(token_map)

    def encode(self, sequence: str, max_codons: int) -> List[int]:
        seq = sequence.upper().replace(" ", "").replace("\n", "")
        codon_count = len(seq) // 3
        if len(seq) % 3 != 0:
            logging.debug("Sequence length not divisible by 3; truncating remainder.")
        codon_count = min(codon_count, max_codons)
        tokens = [self.start_token_id]
        for i in range(codon_count):
            codon = seq[i * 3 : i * 3 + 3]
            token_id = self.token_to_id.get(codon, self.unk_token_id)
            tokens.append(token_id)
        tokens.append(self.end_token_id)
        return tokens


class CodonDataset(Dataset):
    def __init__(self, fasta_path: Path, tokenizer: CodonTokenizer, max_codons: int):
        self.samples: List[List[int]] = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            encoded = tokenizer.encode(str(record.seq), max_codons)
            if len(encoded) > 2:
                self.samples.append(encoded)

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> List[int]:
        return self.samples[idx]


@dataclass
class MaskingConfig:
    mask_prob: float
    random_mask: bool


class CodonMaskingCollator:
    def __init__(self, tokenizer: CodonTokenizer, masking: MaskingConfig, max_length: int):
        self.tokenizer = tokenizer
        self.masking = masking
        self.max_length = max_length

    def __call__(self, batch: List[List[int]]) -> Dict[str, torch.Tensor]:
        batch_ids = []
        for ids in batch:
            if len(ids) > self.max_length:
                ids = ids[: self.max_length - 1] + [self.tokenizer.end_token_id]
            batch_ids.append(torch.tensor(ids, dtype=torch.long))

        input_ids = torch.nn.utils.rnn.pad_sequence(
            batch_ids, batch_first=True, padding_value=self.tokenizer.pad_token_id
        )
        attention_mask = (input_ids != self.tokenizer.pad_token_id).long()
        labels = input_ids.clone()
        labels[labels == self.tokenizer.pad_token_id] = -100

        mask_candidates = (
            (input_ids != self.tokenizer.pad_token_id)
            & (input_ids != self.tokenizer.start_token_id)
            & (input_ids != self.tokenizer.end_token_id)
        )
        mask_positions = mask_candidates

        if mask_positions.any():
            if self.masking.random_mask:
                rand = torch.rand(input_ids.shape)
                mask_token = self.tokenizer.mask_token_id
                random_tokens = torch.tensor(self.tokenizer.codon_ids)
                random_choices = random_tokens[
                    torch.randint(0, len(random_tokens), input_ids.shape)
                ]
                mask_mask = mask_positions & (rand < 0.8)
                random_mask = mask_positions & (rand >= 0.8) & (rand < 0.9)
                input_ids[mask_mask] = mask_token
                input_ids[random_mask] = random_choices[random_mask]
            else:
                input_ids[mask_positions] = self.tokenizer.mask_token_id

        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "labels": labels,
        }


class PerplexityCallback(TrainerCallback):
    def on_evaluate(self, args, state, control, metrics=None, **kwargs):
        if metrics and "eval_loss" in metrics:
            metrics["eval_perplexity"] = math.exp(metrics["eval_loss"])


def load_config(path: Path) -> Dict[str, dict]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def resolve_value(primary, fallback):
    return primary if primary is not None else fallback


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fine-tune BART for codon optimization.")
    parser.add_argument("--pretrained", required=True, help="Path to pretrained checkpoint.")
    parser.add_argument("--config", required=True, help="YAML config path.")
    parser.add_argument("--train", default=None, help="Train FASTA path.")
    parser.add_argument("--val", default=None, help="Validation FASTA path.")
    parser.add_argument("--tokenizer", default=None, help="Codon tokenizer JSON.")
    parser.add_argument("--output", default=None, help="Output directory.")
    parser.add_argument("--steps", type=int, default=None, help="Max training steps.")
    parser.add_argument("--batch-size", type=int, default=None, help="Batch size.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config = load_config(Path(args.config))
    model_cfg = config.get("model", {})
    train_cfg = config.get("training", {})
    mask_cfg = config.get("masking", {})
    data_cfg = config.get("data", {})
    output_cfg = config.get("output", {})
    wandb_cfg = config.get("wandb", {})

    train_path = resolve_value(args.train, data_cfg.get("train"))
    val_path = resolve_value(args.val, data_cfg.get("val"))
    tokenizer_path = resolve_value(args.tokenizer, data_cfg.get("tokenizer"))
    output_dir = resolve_value(args.output, output_cfg.get("save_dir"))

    if not all([train_path, val_path, tokenizer_path, output_dir]):
        logging.error("Missing required paths: train, val, tokenizer, output.")
        return 1

    steps = resolve_value(args.steps, train_cfg.get("steps", 150000))
    batch_size = resolve_value(args.batch_size, train_cfg.get("batch_size", 32))
    warmup_steps = train_cfg.get("warmup_steps", 15000)
    learning_rate = train_cfg.get("learning_rate", 5e-5)
    gradient_accumulation = train_cfg.get("gradient_accumulation_steps", 2)
    fp16 = bool(train_cfg.get("fp16", True))
    weight_decay = train_cfg.get("weight_decay", 0.01)
    max_grad_norm = train_cfg.get("max_grad_norm", 1.0)

    mask_percent = mask_cfg.get("mask_percent", 1.0)
    random_mask = bool(mask_cfg.get("random_mask", False))

    max_codons = data_cfg.get("max_length", 512)
    max_position_embeddings = model_cfg.get("max_position_embeddings", max_codons)
    total_tokens = max_codons + 2
    if max_position_embeddings < total_tokens:
        logging.warning(
            "max_position_embeddings (%d) < max_codons+2 (%d); adjusting upward.",
            max_position_embeddings,
            total_tokens,
        )
        max_position_embeddings = total_tokens

    tokenizer = CodonTokenizer.from_json(Path(tokenizer_path))
    model_config = BartConfig(
        vocab_size=len(tokenizer.token_to_id),
        d_model=model_cfg.get("d_model", 256),
        encoder_layers=model_cfg.get("encoder_layers", 6),
        decoder_layers=model_cfg.get("decoder_layers", 6),
        encoder_attention_heads=model_cfg.get("encoder_attention_heads", 8),
        decoder_attention_heads=model_cfg.get("decoder_attention_heads", 8),
        encoder_ffn_dim=model_cfg.get("encoder_ffn_dim", 256),
        decoder_ffn_dim=model_cfg.get("decoder_ffn_dim", 256),
        max_position_embeddings=max_position_embeddings,
        dropout=model_cfg.get("dropout", 0.1),
        pad_token_id=tokenizer.pad_token_id,
        bos_token_id=tokenizer.start_token_id,
        eos_token_id=tokenizer.end_token_id,
        decoder_start_token_id=tokenizer.start_token_id,
    )

    model = BartForConditionalGeneration.from_pretrained(args.pretrained, config=model_config)

    train_dataset = CodonDataset(Path(train_path), tokenizer, max_codons)
    val_dataset = CodonDataset(Path(val_path), tokenizer, max_codons)
    if len(train_dataset) == 0 or len(val_dataset) == 0:
        logging.error("Training or validation dataset is empty.")
        return 1

    masking = MaskingConfig(mask_prob=mask_percent, random_mask=random_mask)
    collator = CodonMaskingCollator(tokenizer, masking, max_length=total_tokens)

    report_to = ["tensorboard"]
    if wandb_cfg.get("project"):
        try:
            import wandb  # noqa: F401

            report_to.append("wandb")
            os.environ["WANDB_PROJECT"] = str(wandb_cfg.get("project"))
            if wandb_cfg.get("entity"):
                os.environ["WANDB_ENTITY"] = str(wandb_cfg.get("entity"))
        except Exception:
            logging.warning("wandb not available; continuing without WandB logging.")

    training_args = TrainingArguments(
        output_dir=output_dir,
        max_steps=steps,
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        gradient_accumulation_steps=gradient_accumulation,
        evaluation_strategy="steps",
        eval_steps=output_cfg.get("eval_steps", 5000),
        save_strategy="steps",
        save_steps=output_cfg.get("save_steps", 10000),
        logging_steps=output_cfg.get("logging_steps", 100),
        learning_rate=learning_rate,
        warmup_steps=warmup_steps,
        weight_decay=weight_decay,
        max_grad_norm=max_grad_norm,
        fp16=fp16,
        report_to=report_to,
        logging_dir="results/training_logs",
        remove_unused_columns=False,
        run_name=wandb_cfg.get("name"),
        load_best_model_at_end=True,
        metric_for_best_model="eval_loss",
        greater_is_better=False,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        data_collator=collator,
        callbacks=[PerplexityCallback(), EarlyStoppingCallback(early_stopping_patience=5)],
    )

    trainer.train()
    trainer.save_model(output_dir)
    metrics = trainer.evaluate()
    if "eval_loss" in metrics:
        metrics["eval_perplexity"] = math.exp(metrics["eval_loss"])
    trainer.log_metrics("eval", metrics)
    trainer.save_metrics("eval", metrics)
    trainer.save_state()
    return 0


if __name__ == "__main__":
    sys.exit(main())
