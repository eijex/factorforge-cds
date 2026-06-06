"""Security guardrail: raw CDS/protein sequences must not appear in logs or error messages."""
from __future__ import annotations
import logging
import pytest
from factorforge.utils.sequence_validator import validate_cds_output
from factorforge.analysis.metrics import calculate_cai, load_codon_usage_table


LONG_SEQ = "ATGAAACCG" * 30  # 270 bp — long enough to be identifiable as a sequence


def test_validate_cds_output_error_does_not_contain_raw_sequence():
    """Validation errors must not echo back the input sequences."""
    result = validate_cds_output("MMMM", LONG_SEQ)
    if not result["passed"]:
        for err in result.get("errors", []):
            assert LONG_SEQ not in err, f"raw sequence leaked into error: {err[:80]}"


def test_logging_does_not_emit_raw_sequence(caplog):
    """No log record should contain a full raw nucleotide sequence."""
    weights = load_codon_usage_table().codon_weights
    with caplog.at_level(logging.DEBUG, logger="factorforge"):
        calculate_cai(LONG_SEQ, weights)
    for record in caplog.records:
        assert LONG_SEQ not in record.getMessage(), (
            f"raw sequence leaked into log ({record.name}): {record.getMessage()[:80]}"
        )


def test_exception_traceback_does_not_leak_sequence():
    """Custom exceptions raised by factorforge must not embed raw sequences in message or traceback.
    Python traceback can dump function arguments — this guards against that path."""
    import traceback as tb_mod
    try:
        # Mismatched protein length triggers validation failure; captures any exception path
        validate_cds_output("M" * 90, LONG_SEQ)
    except Exception as exc:
        tb_str = tb_mod.format_exc()
        assert LONG_SEQ not in str(exc), "raw sequence leaked in exception message"
        assert LONG_SEQ not in tb_str, "raw sequence leaked in traceback dump"
    # If function returns a result dict (no exception), check error strings in result
    result = validate_cds_output("M" * 90, LONG_SEQ)
    for err in result.get("errors", []):
        assert LONG_SEQ not in err, f"raw sequence leaked in validation error: {err[:80]}"
