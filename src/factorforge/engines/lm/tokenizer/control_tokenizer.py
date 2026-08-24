"""FactorForge Control Tokenizer (v3.5.0 ML Track).

Defines the control tokens, special tokens, and codon vocabulary for
FactorForge-LM (SynCodonLM-V2), enabling host-aware and constraint-guided
codon sequence encoding and decoding.
"""

from typing import List, Dict

# Special Tokens
PAD_TOKEN = "<pad>"
UNK_TOKEN = "<unk>"
BOS_TOKEN = "<s>"
EOS_TOKEN = "</s>"
MASK_TOKEN = "<mask>"

SPECIAL_TOKENS = [PAD_TOKEN, UNK_TOKEN, BOS_TOKEN, EOS_TOKEN, MASK_TOKEN]

# Host Control Tokens
HOST_TOKENS = [
    "<host:nbenthamiana>",
    "<host:ntabacum>",
    "<host:by2>",
]

# GC Band Control Tokens
GC_TOKENS = [
    "<gc:40-47>",
    "<gc:47-55>",
    "<gc:55-65>",
]

# Assembly Control Tokens
ASSEMBLY_TOKENS = [
    "<type2is:clean>",
    "<moclo:level0>",
]

# Standard Genetic Code Codons (61 Sense + 3 Stop)
SENSE_CODONS = [
    "GCT",
    "GCC",
    "GCA",
    "GCG",  # Ala (A)
    "CGT",
    "CGC",
    "CGA",
    "CGG",
    "AGA",
    "AGG",  # Arg (R)
    "AAT",
    "AAC",  # Asn (N)
    "GAT",
    "GAC",  # Asp (D)
    "TGT",
    "TGC",  # Cys (C)
    "CAA",
    "CAG",  # Gln (Q)
    "GAA",
    "GAG",  # Glu (E)
    "GGT",
    "GGC",
    "GGA",
    "GGG",  # Gly (G)
    "CAT",
    "CAC",  # His (H)
    "ATT",
    "ATC",
    "ATA",  # Ile (I)
    "TTA",
    "TTG",
    "CTA",
    "CTC",
    "CTG",
    "CTT",  # Leu (L)
    "AAA",
    "AAG",  # Lys (K)
    "ATG",  # Met (M)
    "TTT",
    "TTC",  # Phe (F)
    "CCT",
    "CCC",
    "CCA",
    "CCG",  # Pro (P)
    "TCT",
    "TCC",
    "TCA",
    "TCG",
    "AGT",
    "AGC",  # Ser (S)
    "ACT",
    "ACC",
    "ACA",
    "ACG",  # Thr (T)
    "TGG",  # Trp (W)
    "TAT",
    "TAC",  # Tyr (Y)
    "GTT",
    "GTC",
    "GTA",
    "GTG",  # Val (V)
]

STOP_CODONS = ["TAA", "TAG", "TGA"]
ALL_CODONS = SENSE_CODONS + STOP_CODONS

# Amino Acid Tokens for Encoder Input
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY*")


class FactorForgeControlTokenizer:
    """Tokenizer for FactorForge-LM supporting control tokens and codons."""

    def __init__(self):
        self.vocab: List[str] = []
        self.token_to_id: Dict[str, int] = {}
        self.id_to_token: Dict[int, str] = {}
        self._build_vocab()

    def _build_vocab(self):
        """Constructs vocabulary mappings."""
        all_tokens = (
            SPECIAL_TOKENS + HOST_TOKENS + GC_TOKENS + ASSEMBLY_TOKENS + AMINO_ACIDS + ALL_CODONS
        )
        seen = set()
        for tok in all_tokens:
            if tok not in seen:
                seen.add(tok)
                self.vocab.append(tok)

        self.token_to_id = {tok: idx for idx, tok in enumerate(self.vocab)}
        self.id_to_token = {idx: tok for idx, tok in enumerate(self.vocab)}

    @property
    def vocab_size(self) -> int:
        return len(self.vocab)

    @property
    def pad_id(self) -> int:
        return self.token_to_id[PAD_TOKEN]

    @property
    def unk_id(self) -> int:
        return self.token_to_id[UNK_TOKEN]

    @property
    def bos_id(self) -> int:
        return self.token_to_id[BOS_TOKEN]

    @property
    def eos_id(self) -> int:
        return self.token_to_id[EOS_TOKEN]

    @property
    def mask_id(self) -> int:
        return self.token_to_id[MASK_TOKEN]

    def encode_encoder_input(
        self,
        amino_acids: str,
        host: str = "nbenthamiana",
        gc_band: str = "40-47",
        type2is_clean: bool = True,
    ) -> List[int]:
        """Encodes an amino acid sequence with prefix control tokens for the Encoder."""
        tokens = [BOS_TOKEN]

        # Add Host Token
        host_tok = f"<host:{host.lower()}>"
        tokens.append(host_tok if host_tok in self.token_to_id else "<host:nbenthamiana>")

        # Add GC Token
        gc_tok = f"<gc:{gc_band}>"
        tokens.append(gc_tok if gc_tok in self.token_to_id else "<gc:40-47>")

        # Add Assembly Token
        if type2is_clean:
            tokens.append("<type2is:clean>")

        # Add Amino Acid Tokens
        for aa in amino_acids.upper():
            tokens.append(aa if aa in self.token_to_id else UNK_TOKEN)

        tokens.append(EOS_TOKEN)
        return [self.token_to_id.get(t, self.unk_id) for t in tokens]

    def decode_ids(self, ids: List[int], skip_special_tokens: bool = True) -> str:
        """Decodes token IDs back into a sequence string."""
        tokens = []
        for idx in ids:
            tok = self.id_to_token.get(idx, UNK_TOKEN)
            if skip_special_tokens and tok in SPECIAL_TOKENS:
                continue
            tokens.append(tok)
        return "".join(tokens)


if __name__ == "__main__":
    tokenizer = FactorForgeControlTokenizer()
    print(f"FactorForge Control Tokenizer initialized. Vocab size: {tokenizer.vocab_size}")
    encoded = tokenizer.encode_encoder_input("MAKW", host="nbenthamiana", gc_band="40-47")
    print(f"Sample Encoded Input IDs: {encoded}")
    print(f"Sample Decoded Tokens: {[tokenizer.id_to_token[i] for i in encoded]}")
