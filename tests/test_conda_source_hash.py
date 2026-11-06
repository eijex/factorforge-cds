"""Published-source binding cannot silently leave a pending recipe hash."""

import importlib.util
from pathlib import Path

import pytest

path = Path(__file__).resolve().parents[1] / "scripts/update_conda_recipe_sha.py"
spec = importlib.util.spec_from_file_location("conda_source", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.mark.parametrize("old", ["REPLACE_AFTER_PYPI_3_5_12_PUBLICATION", "a" * 64])
def test_bind_real_digest_over_placeholder_or_old_hash(old):
    text = '{% set version = "3.5.11" %}\nsource:\n  sha256: ' + old + '\n'
    result = module.replace_recipe_source(text, "3.5.12", "b" * 64)
    assert 'version = "3.5.12"' in result
    assert "sha256: " + "b" * 64 in result
    assert old not in result


@pytest.mark.parametrize("text", ["", "sha256: x\nsha256: y\n"])
def test_reject_missing_or_ambiguous_recipe(text):
    with pytest.raises(ValueError):
        module.replace_recipe_source(text, "3.5.12", "b" * 64)


def test_reject_invalid_digest():
    with pytest.raises(ValueError):
        module.replace_recipe_source("", "3.5.12", "not-a-digest")
