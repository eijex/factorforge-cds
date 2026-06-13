"""Tests for the sfGFP worked example and ValidationHub public schema.

Covers:
- Input sequence integrity (AC1, AC2)
- Frozen output consistency — design_package.json (AC3, AC4)
- Frozen output consistency — validation_summary.json (AC4, AC5, AC6)
- Scoring contract v1.1 semantics (AC5, AC6)
- Public schema validity (AC7, AC8)
- validation_summary validates against public schema (AC8)
- Run example reproducibility (AC2)
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIR = ROOT / "examples" / "worked_example"
OUTPUT_DIR = EXAMPLE_DIR / "output"
SCHEMA_DIR = ROOT / "docs" / "validationhub"

GC_MIN = 55.0
GC_MAX = 65.0


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# TC1 — sfGFP input sequence
# ---------------------------------------------------------------------------

def test_sfgfp_sequence_header():
    faa = EXAMPLE_DIR / "input_sequence.faa"
    first_line = faa.read_text(encoding="utf-8").strip().split("\n")[0]
    assert first_line.startswith(">sfGFP_worked_example"), "FASTA header mismatch"
    assert "Pédelacq" in first_line or "Pedelacq" in first_line, "citation missing"
    assert "PDB:2B3P" in first_line, "PDB source missing"


def test_sfgfp_sequence_length():
    faa = EXAMPLE_DIR / "input_sequence.faa"
    lines = faa.read_text(encoding="utf-8").strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    assert len(seq) == 236, f"expected 236 aa, got {len(seq)}"


def test_sfgfp_sequence_no_x_residues():
    faa = EXAMPLE_DIR / "input_sequence.faa"
    lines = faa.read_text(encoding="utf-8").strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    assert "X" not in seq, "X (ambiguous) residue must not appear in worked example"


# ---------------------------------------------------------------------------
# TC2 — Reproducibility
# ---------------------------------------------------------------------------

def test_run_example_deterministic():
    """run_example.py (no --freeze) must exit 0 when frozen outputs exist."""
    result = subprocess.run(
        [sys.executable, str(EXAMPLE_DIR / "run_example.py")],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )
    assert result.returncode == 0, (
        f"run_example.py exited {result.returncode}\nstdout:{result.stdout}\nstderr:{result.stderr}"
    )
    assert "OK" in result.stdout, "expected 'OK' in output"


# ---------------------------------------------------------------------------
# TC3 — design_package.json claim boundary
# ---------------------------------------------------------------------------

def test_design_package_claim_boundary():
    dp = _load(OUTPUT_DIR / "design_package.json")
    cb = dp["claim_boundary"]
    assert cb["in_silico_only"] is True
    assert cb["no_yield_claim"] is True
    assert cb["no_wet_lab_claim"] is True
    assert cb["no_clinical_claim"] is True


def test_design_package_profile():
    dp = _load(OUTPUT_DIR / "design_package.json")
    assert dp["optimization"]["profile"] == "assembly_friendly"
    assert dp["optimization"]["engine"] == "profile"


# ---------------------------------------------------------------------------
# TC4 — validation_summary.json scoring contract version
# ---------------------------------------------------------------------------

def test_validation_summary_scoring_contract():
    vs = _load(OUTPUT_DIR / "validation_summary.json")
    assert vs["computational"]["scoring_contract_version"] == "v1.1", (
        "scoring_contract_version must be v1.1"
    )


def test_validation_summary_definition_const():
    vs = _load(OUTPUT_DIR / "validation_summary.json")
    defn = vs["computational"]["multi_constraint_pass_definition"]
    assert defn == "biological_pass AND assembly_pass AND gc_in_target_range", (
        f"unexpected definition: {defn}"
    )


# ---------------------------------------------------------------------------
# TC5 — gc_in_target_range is primitive boolean
# ---------------------------------------------------------------------------

def test_gc_in_target_range_is_bool():
    vs = _load(OUTPUT_DIR / "validation_summary.json")
    assert isinstance(vs["computational"]["gc_in_target_range"], bool)


def test_gc_in_target_range_consistent_with_design_package():
    dp = _load(OUTPUT_DIR / "design_package.json")
    vs = _load(OUTPUT_DIR / "validation_summary.json")
    gc = dp["metrics"]["gc_percent"]
    expected = GC_MIN <= gc <= GC_MAX
    assert vs["computational"]["gc_in_target_range"] == expected, (
        f"gc_in_target_range mismatch: gc={gc}, expected={expected}"
    )


# ---------------------------------------------------------------------------
# TC6 — multi_constraint_pass is derived from primitives
# ---------------------------------------------------------------------------

def test_multi_constraint_pass_is_derived():
    dp = _load(OUTPUT_DIR / "design_package.json")
    vs = _load(OUTPUT_DIR / "validation_summary.json")
    gc = dp["metrics"]["gc_percent"]
    bio = dp["validation"]["biological_pass"]
    asm = dp["validation"]["assembly_pass"]
    gc_ok = GC_MIN <= gc <= GC_MAX
    expected = bio and asm and gc_ok
    assert vs["computational"]["multi_constraint_pass"] == expected, (
        f"multi_constraint_pass ({vs['computational']['multi_constraint_pass']}) "
        f"!= derived ({expected}); bio={bio}, asm={asm}, gc_ok={gc_ok}"
    )


# ---------------------------------------------------------------------------
# TC7 — public schema is valid JSON Schema (Draft 2020-12)
# ---------------------------------------------------------------------------

def test_public_schema_is_valid_json_schema():
    jsonschema = pytest.importorskip("jsonschema")
    schema = _load(SCHEMA_DIR / "intake_schema_v0.1.public.json")
    jsonschema.Draft202012Validator.check_schema(schema)


# ---------------------------------------------------------------------------
# TC8 — validation_summary.json validates against public schema
# ---------------------------------------------------------------------------

def test_validation_summary_validates_against_public_schema():
    jsonschema = pytest.importorskip("jsonschema")
    schema = _load(SCHEMA_DIR / "intake_schema_v0.1.public.json")
    vs = _load(OUTPUT_DIR / "validation_summary.json")
    validator = jsonschema.Draft202012Validator(schema)
    errors = list(validator.iter_errors(vs))
    assert not errors, f"Schema validation failed:\n" + "\n".join(str(e) for e in errors)
