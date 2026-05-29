"""
Optional ML-based scoring helpers for the profile engine.

SynCodonLM is loaded lazily so the default FactorForge installation remains
free of torch/transformers runtime requirements.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

SYNCODONLM_MODEL_NAME = "jheuschkel/SynCodonLM-V2"
DEFAULT_SPECIES_TOKEN_TYPE = 500
NEUTRAL_SYNCODONLM_SCORE = 0.5

# SynCodonLM's README lists model-organism token type IDs. N. benthamiana is
# not listed there, so it intentionally falls back to the model default.
SPECIES_TOKEN_TYPES = {
    "e_coli": 67,
    "escherichia_coli": 67,
    "s_cerevisiae": 108,
    "saccharomyces_cerevisiae": 108,
    "c_elegans": 187,
    "caenorhabditis_elegans": 187,
    "d_melanogaster": 178,
    "drosophila_melanogaster": 178,
    "d_rerio": 468,
    "danio_rerio": 468,
    "m_musculus": 321,
    "mus_musculus": 321,
    "a_thaliana": 266,
    "arabidopsis_thaliana": 266,
    "h_sapiens": 317,
    "homo_sapiens": 317,
    "c_griseus": 394,
    "cricetulus_griseus": 394,
}

_DEFAULT_SCORERS: dict[tuple[str, str | None], SynCodonLMScorer] = {}


def _organism_key(organism: str) -> str:
    """Normalize organism names for token type lookup."""
    return organism.strip().lower().replace(".", "").replace(" ", "_").replace("-", "_")


def _species_token_type(organism: str) -> int:
    """Return the SynCodonLM species token type or the model default."""
    return SPECIES_TOKEN_TYPES.get(_organism_key(organism), DEFAULT_SPECIES_TOKEN_TYPE)


def _split_codons(sequence: str) -> list[str]:
    """Clean a CDS sequence and return codon tokens."""
    cleaned = sequence.upper().replace("U", "T").replace(" ", "").replace("\n", "")
    cleaned = cleaned.replace("\r", "").replace("\t", "")
    if not cleaned:
        raise ValueError("sequence must not be empty")
    invalid = sorted({base for base in cleaned if base not in {"A", "T", "G", "C"}})
    if invalid:
        raise ValueError(f"invalid nucleotide(s): {', '.join(invalid)}")
    if len(cleaned) % 3 != 0:
        raise ValueError("sequence length must be divisible by 3")
    return [cleaned[i : i + 3] for i in range(0, len(cleaned), 3)]


@dataclass
class SynCodonLMScorer:
    """Lazy SynCodonLM masked-language-model scorer."""

    model_name: str = SYNCODONLM_MODEL_NAME
    device: str | None = None
    _tokenizer: Any = field(init=False, default=None, repr=False)
    _model: Any = field(init=False, default=None, repr=False)
    _torch: Any = field(init=False, default=None, repr=False)
    _mask_id: int | None = field(init=False, default=None, repr=False)
    _bos_id: int | None = field(init=False, default=None, repr=False)
    _eos_id: int | None = field(init=False, default=None, repr=False)
    _supports_token_type_ids: bool = field(init=False, default=False, repr=False)
    _available: bool | None = field(init=False, default=None, repr=False)

    def score(self, sequence: str, organism: str = "Nicotiana_benthamiana") -> float:
        """Score a CDS with SynCodonLM.

        Returns a 0-1 geometric mean probability for the observed codons under
        one-codon-at-a-time masking. If the optional ML stack is unavailable or
        scoring fails, returns a neutral 0.5.
        """
        if not self._load():
            return NEUTRAL_SYNCODONLM_SCORE

        try:
            codons = _split_codons(sequence)
            token_ids = self._codon_token_ids(codons)
            if not token_ids:
                return NEUTRAL_SYNCODONLM_SCORE

            log_probs: list[float] = []
            with self._torch.no_grad():
                for index, target_id in enumerate(token_ids):
                    inputs = self._masked_inputs(token_ids, index, _species_token_type(organism))
                    logits = self._model(**inputs).logits.squeeze(0)
                    row = logits[index + 1]
                    log_prob = self._torch.log_softmax(row, dim=-1)[target_id]
                    log_probs.append(float(log_prob.item()))

            mean_log_prob = sum(log_probs) / len(log_probs)
            return round(max(0.0, min(1.0, math.exp(mean_log_prob))), 3)
        except Exception as exc:
            logger.warning("SynCodonLM scoring failed; returning neutral score: %s", exc)
            return NEUTRAL_SYNCODONLM_SCORE

    def _load(self) -> bool:
        """Load torch, tokenizer, and model on first use."""
        if self._available is not None:
            return self._available

        try:
            import torch
            from transformers import AutoConfig, AutoModelForMaskedLM, AutoTokenizer
        except ImportError as exc:
            logger.warning(
                "SynCodonLM scoring requires optional dependencies; install factorforge-cds[ml]. "
                "Returning neutral score. Missing dependency: %s",
                exc,
            )
            self._available = False
            return False

        try:
            device = torch.device(self.device) if self.device else torch.device("cpu")
            tokenizer = AutoTokenizer.from_pretrained(self.model_name, use_fast=True)
            config = AutoConfig.from_pretrained(self.model_name)
            model = AutoModelForMaskedLM.from_pretrained(self.model_name, config=config)
            model.to(device).eval()

            mask_id = tokenizer.mask_token_id
            bos_id = tokenizer.bos_token_id or tokenizer.cls_token_id
            eos_id = tokenizer.eos_token_id or tokenizer.sep_token_id
            if mask_id is None or bos_id is None or eos_id is None:
                raise ValueError("SynCodonLM tokenizer requires BOS/EOS and MASK tokens")

            self._torch = torch
            self._tokenizer = tokenizer
            self._model = model
            self._mask_id = int(mask_id)
            self._bos_id = int(bos_id)
            self._eos_id = int(eos_id)
            self._supports_token_type_ids = (
                "token_type_ids" in getattr(tokenizer, "model_input_names", [])
                and "NoTokenType" not in getattr(config, "name_or_path", "")
            )
            self._available = True
            return True
        except Exception as exc:
            logger.warning("Could not load SynCodonLM model; returning neutral score: %s", exc)
            self._available = False
            return False

    def _codon_token_ids(self, codons: list[str]) -> list[int]:
        """Convert codon strings into tokenizer IDs."""
        unk_id = self._tokenizer.unk_token_id
        token_ids: list[int] = []
        for codon in codons:
            token_id = self._tokenizer.convert_tokens_to_ids(codon)
            if token_id is None or token_id == unk_id:
                raise ValueError(f"codon token is not in SynCodonLM vocabulary: {codon}")
            token_ids.append(int(token_id))
        return token_ids

    def _masked_inputs(
        self,
        token_ids: list[int],
        mask_index: int,
        species_token_type: int,
    ) -> dict[str, Any]:
        """Build a single-position masked input tensor."""
        body = list(token_ids)
        body[mask_index] = self._mask_id
        input_ids = self._torch.tensor(
            [[self._bos_id] + body + [self._eos_id]],
            dtype=self._torch.long,
            device=self._model.device,
        )
        inputs = {
            "input_ids": input_ids,
            "attention_mask": self._torch.ones_like(input_ids),
        }
        if self._supports_token_type_ids:
            inputs["token_type_ids"] = self._torch.full_like(input_ids, int(species_token_type))
        return inputs


def calculate_syncodonlm_score(
    sequence: str,
    organism: str = "Nicotiana_benthamiana",
) -> float:
    """Calculate an optional SynCodonLM score for a CDS sequence."""
    scorer_key = (SYNCODONLM_MODEL_NAME, None)
    scorer = _DEFAULT_SCORERS.get(scorer_key)
    if scorer is None:
        scorer = SynCodonLMScorer()
        _DEFAULT_SCORERS[scorer_key] = scorer
    return scorer.score(sequence, organism=organism)
