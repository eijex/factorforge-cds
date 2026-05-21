"""Tests for Run 4 v2 pseudo-label generation."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = ROOT / "scripts" / "1_data_preparation"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from generate_v2_pseudolabels import generate_v2_pseudolabels, split_jsonl  # noqa: E402


def test_generate_v2_pseudolabels_writes_valid_jsonl(tmp_path: Path) -> None:
    fasta = tmp_path / "proteins.fasta"
    output = tmp_path / "labels.jsonl"
    fasta.write_text(">p1\nMAF\n>p2\nMST\n", encoding="utf-8")

    counts = generate_v2_pseudolabels(fasta, output, profile="high_cai")
    rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]

    assert counts["read"] == 2
    assert counts["written"] == 2
    assert rows[0]["protein_id"] == "p1"
    assert rows[0]["sequence"] == "MAF"
    assert rows[0]["protein_sequence"] == "MAF"
    assert rows[0]["amino_acid_sequence"] == "MAF"
    assert rows[0]["dna_sequence"] == rows[0]["dna_seq"]
    assert rows[0]["validator"]["passed"] is True
    assert len(rows[0]["dna_seq"]) == len(rows[0]["sequence"]) * 3
    assert len(rows[0]["codon_sequence"].split()) == len(rows[0]["sequence"])


def test_generate_v2_pseudolabels_accepts_cds_fasta_and_splits(tmp_path: Path) -> None:
    fasta = tmp_path / "cds.fasta"
    output = tmp_path / "labels.jsonl"
    train = tmp_path / "train.jsonl"
    eval_file = tmp_path / "eval.jsonl"
    fasta.write_text(">p1\nATGGCTTTCTAA\n>p2\nATGTCTACTTGA\n", encoding="utf-8")

    counts = generate_v2_pseudolabels(fasta, output, profile="high_cai", input_kind="cds")
    split = split_jsonl(output, train, eval_file, train_split=0.5, seed=1)
    rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]

    assert counts["read"] == 2
    assert counts["written"] == 2
    assert rows[0]["amino_acid_sequence"] == "MAF"
    assert rows[1]["amino_acid_sequence"] == "MST"
    assert split == {"total": 2, "train": 1, "eval": 1}
