"""
Machine-learning-based codon optimization for plants.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import torch
from transformers import BartForConditionalGeneration

CODON_TO_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


class CodonTokenizer:
    def __init__(self, token_map: Dict[str, int]):
        self.token_to_id = token_map
        self.id_to_token = {idx: token for token, idx in token_map.items()}
        self.pad_token_id = token_map["[PAD]"]
        self.unk_token_id = token_map["[UNK]"]
        self.mask_token_id = token_map["[MASK]"]
        self.start_token_id = token_map["[START]"]
        self.end_token_id = token_map["[END]"]

    @classmethod
    def from_json(cls, path: Path) -> "CodonTokenizer":
        with path.open("r", encoding="utf-8") as handle:
            token_map = json.load(handle)
        required = ["[PAD]", "[UNK]", "[MASK]", "[START]", "[END]"]
        missing = [token for token in required if token not in token_map]
        if missing:
            raise ValueError(f"Tokenizer missing special tokens: {missing}")
        return cls(token_map)

    def encode_dna(self, dna_seq: str) -> List[int]:
        seq = dna_seq.upper()
        tokens = [self.start_token_id]
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if len(codon) != 3:
                continue
            token_id = self.token_to_id.get(codon, self.unk_token_id)
            tokens.append(token_id)
        tokens.append(self.end_token_id)
        return tokens

    def decode_dna(self, ids: List[int]) -> str:
        codons: List[str] = []
        for idx in ids:
            token = self.id_to_token.get(int(idx))
            if token and len(token) == 3 and all(base in "ACGT" for base in token):
                codons.append(token)
        return "".join(codons)


class PlantCodonOptimizer:
    """
    Optimize codon usage using a trained BART model and codon frequency table.
    """

    def __init__(
        self,
        model_path: str,
        codon_table_path: str,
        tokenizer_path: str,
        organism: str = "N.benthamiana",
    ) -> None:
        self.organism = organism
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = BartForConditionalGeneration.from_pretrained(model_path).to(self.device)
        self.model.eval()

        self.codon_table = self._load_codon_table(Path(codon_table_path))
        self.tokenizer = CodonTokenizer.from_json(Path(tokenizer_path))
        self.codon_weights = self._build_codon_weights(self.codon_table)
        self.best_codon_for_aa = self._best_codon_map(self.codon_table)

    def optimize(self, protein_sequence: str, beam_size: int = 5) -> str:
        """
        Generate an optimized DNA sequence for a protein input.
        """
        protein_seq = self._normalize_protein(protein_sequence)
        baseline_dna = self._reverse_translate(protein_seq)
        input_ids = torch.tensor(
            [self.tokenizer.encode_dna(baseline_dna)],
            dtype=torch.long,
            device=self.device,
        )
        attention_mask = (input_ids != self.tokenizer.pad_token_id).long()
        max_length = input_ids.shape[1]

        with torch.no_grad():
            outputs = self.model.generate(
                input_ids=input_ids,
                attention_mask=attention_mask,
                num_beams=beam_size,
                max_length=max_length,
                min_length=max_length,
                early_stopping=True,
            )
        return self.tokenizer.decode_dna(outputs[0].tolist())

    def calculate_cai(self, dna_sequence: str) -> float:
        """
        Calculate Codon Adaptation Index using codon frequencies.
        """
        seq = dna_sequence.upper()
        codon_count = len(seq) // 3
        if codon_count == 0:
            return 0.0
        weights = []
        for i in range(codon_count):
            codon = seq[i * 3 : i * 3 + 3]
            weight = self.codon_weights.get(codon, 0.0)
            if weight <= 0:
                return 0.0
            weights.append(weight)
        log_sum = sum(math.log(w) for w in weights)
        return math.exp(log_sum / len(weights))

    def compare_sequences(self, original_dna: str, optimized_dna: str) -> Dict[str, float]:
        """
        Compare codon-by-codon and return metrics.
        """
        original = original_dna.upper()
        optimized = optimized_dna.upper()
        total_codons = min(len(original), len(optimized)) // 3
        changed = 0
        for i in range(total_codons):
            o_codon = original[i * 3 : i * 3 + 3]
            n_codon = optimized[i * 3 : i * 3 + 3]
            if o_codon != n_codon:
                changed += 1

        original_cai = self.calculate_cai(original)
        optimized_cai = self.calculate_cai(optimized)
        return {
            "total_codons": total_codons,
            "changed_codons": changed,
            "change_rate": (changed / total_codons * 100) if total_codons else 0.0,
            "original_cai": original_cai,
            "optimized_cai": optimized_cai,
            "cai_improvement": optimized_cai - original_cai,
            "original_gc": self._gc_content(original),
            "optimized_gc": self._gc_content(optimized),
        }

    def generate_report(
        self,
        protein_name: str,
        protein_seq: str,
        original_dna: str,
        optimized_dna: str,
        output_path: str,
    ) -> None:
        """
        Generate a formatted text report for the optimization.
        """
        metrics = self.compare_sequences(original_dna, optimized_dna)
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8") as handle:
            handle.write("Codon Optimization Report\n")
            handle.write("=" * 60 + "\n")
            handle.write(f"Organism: {self.organism}\n")
            handle.write(f"Protein: {protein_name}\n")
            handle.write(f"Protein length: {len(protein_seq)} aa\n\n")
            handle.write("Metrics\n")
            handle.write("-" * 60 + "\n")
            handle.write(f"Total codons: {metrics['total_codons']}\n")
            handle.write(f"Changed codons: {metrics['changed_codons']}\n")
            handle.write(f"Change rate: {metrics['change_rate']:.2f}%\n")
            handle.write(f"Original CAI: {metrics['original_cai']:.4f}\n")
            handle.write(f"Optimized CAI: {metrics['optimized_cai']:.4f}\n")
            handle.write(f"CAI improvement: {metrics['cai_improvement']:.4f}\n")
            handle.write(f"Original GC: {metrics['original_gc']:.2f}%\n")
            handle.write(f"Optimized GC: {metrics['optimized_gc']:.2f}%\n\n")
            handle.write("Original DNA\n")
            handle.write("-" * 60 + "\n")
            handle.write(original_dna + "\n\n")
            handle.write("Optimized DNA\n")
            handle.write("-" * 60 + "\n")
            handle.write(optimized_dna + "\n")

    def _load_codon_table(self, path: Path) -> Dict[str, float]:
        if not path.exists():
            raise FileNotFoundError(f"Codon table not found: {path}")
        df = pd.read_csv(path)
        columns = {col.lower(): col for col in df.columns}
        if "codon" not in columns or "frequency" not in columns:
            raise ValueError("Codon table must have columns: Codon, Frequency")
        codon_col = columns["codon"]
        freq_col = columns["frequency"]
        codon_freq: Dict[str, float] = {}
        for _, row in df.iterrows():
            codon = str(row[codon_col]).strip().upper()
            try:
                freq = float(row[freq_col])
            except (TypeError, ValueError):
                freq = 0.0
            if len(codon) == 3:
                codon_freq[codon] = freq
        return codon_freq

    def _build_codon_weights(self, codon_freq: Dict[str, float]) -> Dict[str, float]:
        by_aa: Dict[str, List[Tuple[str, float]]] = {}
        for codon, freq in codon_freq.items():
            aa = CODON_TO_AA.get(codon, "*")
            if aa == "*":
                continue
            by_aa.setdefault(aa, []).append((codon, freq))

        weights: Dict[str, float] = {}
        for aa, codons in by_aa.items():
            max_freq = max(freq for _, freq in codons) if codons else 0.0
            for codon, freq in codons:
                weights[codon] = freq / max_freq if max_freq > 0 else 0.0
        return weights

    def _best_codon_map(self, codon_freq: Dict[str, float]) -> Dict[str, str]:
        best: Dict[str, Tuple[str, float]] = {}
        for codon, freq in codon_freq.items():
            aa = CODON_TO_AA.get(codon, "*")
            if aa == "*":
                continue
            current = best.get(aa)
            if current is None or freq > current[1]:
                best[aa] = (codon, freq)
        return {aa: codon for aa, (codon, _) in best.items()}

    def _reverse_translate(self, protein_seq: str) -> str:
        codons = []
        for aa in protein_seq:
            codon = self.best_codon_for_aa.get(aa)
            if codon is None:
                raise ValueError(f"No codon mapping for amino acid: {aa}")
            codons.append(codon)
        return "".join(codons)

    def _normalize_protein(self, protein_sequence: str) -> str:
        seq = protein_sequence.strip().replace("\n", "").replace(" ", "").upper()
        valid = set("ACDEFGHIKLMNPQRSTVWY")
        if not seq:
            raise ValueError("Protein sequence is empty.")
        invalid = {ch for ch in seq if ch not in valid}
        if invalid:
            raise ValueError(f"Invalid amino acids found: {''.join(sorted(invalid))}")
        return seq

    @staticmethod
    def _gc_content(dna_sequence: str) -> float:
        seq = dna_sequence.upper()
        if not seq:
            return 0.0
        gc = seq.count("G") + seq.count("C")
        return (gc / len(seq)) * 100.0
