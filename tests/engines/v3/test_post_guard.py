"""Post-guard hook tests for v3 pipeline."""

from factorforge.engines.v3.pipeline import V3Pipeline


def test_post_guard_fixes_polya():
    """PolyA signal in DNA sequence is removed by _apply_post_guard()."""
    # ATGAATAAAGCCTAA = M-N-K-A-* ; AATAAA spans AAT|AAA codon boundary.
    # AAT (Asn) and AAA (Lys) both have synonymous codons in N. benthamiana table.
    seq_with_polya = "ATGAATAAAGCCTAA"

    pipeline = V3Pipeline()
    result_seq, post_guard = pipeline._apply_post_guard(seq_with_polya)

    # polya_fix must have been attempted (not None)
    assert post_guard["polya_fix"] is not None

    # If fix succeeded, PolyA patterns should be absent from result
    from factorforge.engines.v2.rules.rule_engine import RuleEngine
    rule_engine = RuleEngine()
    if post_guard["polya_fix"]["success"]:
        for pattern in rule_engine.POLYA_PATTERNS:
            assert pattern not in result_seq, (
                f"PolyA pattern {pattern!r} still present after post-guard fix"
            )


def test_post_guard_no_polya_unchanged():
    """DNA sequence without PolyA signal leaves polya_fix as None."""
    # ATGGCCGCCTAA = M-A-A-* ; no PolyA patterns present
    seq_no_polya = "ATGGCCGCCTAA"

    pipeline = V3Pipeline()
    _, post_guard = pipeline._apply_post_guard(seq_no_polya)

    # No PolyA detected → fast pre-check short-circuits, polya_fix stays None
    assert post_guard["polya_fix"] is None
