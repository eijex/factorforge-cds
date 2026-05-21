"""Tokenizer contract tests for v3."""

from factorforge.engines.v3.tokenizer import (
    AA_TOKENS,
    CODON_TOKENS,
    SPECIAL_TOKENS,
    AATokenizer,
    CodonTokenizer,
)


def test_aa_tokenizer_roundtrip():
    tokenizer = AATokenizer.default()
    sequence = "".join(AA_TOKENS)
    encoded = tokenizer.encode(sequence, add_special_tokens=False)
    decoded = tokenizer.decode(encoded, skip_special_tokens=True)
    assert decoded == sequence


def test_codon_tokenizer_roundtrip():
    tokenizer = CodonTokenizer.default()
    sequence = "ATGTTTCCCGGG"
    encoded = tokenizer.encode(sequence, add_special_tokens=False)
    decoded = tokenizer.decode(encoded, skip_special_tokens=True)
    assert decoded == sequence


def test_tokenizer_mapping_persistence(tmp_path):
    tokenizer = CodonTokenizer.default()
    path = tmp_path / "codon_tokenizer.json"
    tokenizer.save(path)
    loaded = CodonTokenizer.load(path)

    assert loaded.token_to_id == tokenizer.token_to_id
    assert loaded.pad_token_id == tokenizer.pad_token_id
    assert loaded.bos_token_id == tokenizer.bos_token_id
    assert loaded.eos_token_id == tokenizer.eos_token_id


def test_special_token_ids_are_stable():
    tokenizer = AATokenizer.default()
    assert tokenizer.token_to_id[SPECIAL_TOKENS[0]] == 0
    assert tokenizer.token_to_id[SPECIAL_TOKENS[1]] == 1
    assert tokenizer.token_to_id[SPECIAL_TOKENS[2]] == 2
    assert tokenizer.token_to_id[SPECIAL_TOKENS[3]] == 3
    assert tokenizer.token_to_id[SPECIAL_TOKENS[4]] == 4
    assert len(CODON_TOKENS) == 64
