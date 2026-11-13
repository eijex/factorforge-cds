"""Tests for FactorForge Default SOP refinement (5' SP lock + Arg/Ser balancing)."""

from __future__ import annotations

from Bio.Seq import Seq

from factorforge.codon.sop_balancer import (
    ARG_CODONS,
    SER_CODONS,
    apply_default_sop_refinement,
)


def test_preserve_5p_signal_peptide_and_balance_arg_ser() -> None:
    # 20-aa SP with non-default synonymous codons + mature region with skewed Arg/Ser and tandem repeats
    orig_sp = "ATGGCCAAAACCAATCTTTTTCTTTTTTTGATTTTTTCTCTTCTTCTTTCTCTTTCTTCTGCT"[:60]
    assert len(orig_sp) == 60
    # Mature: 6 Arg (all AGA, including tandem pair) + 10 Ser (all TCT, including tandem pairs) + Ala spacers
    mature_skewed = (
        ("AGA" * 2) + "GCT" + "AGA" + "GCT" + "AGA" + "GCT" + ("AGA" * 2)
        + "GCT" + ("TCT" * 2) + "GCT" + ("TCT" * 2) + "GCT" + ("TCT" * 2) + "GCT" + ("TCT" * 2) + "GCT" + ("TCT" * 2)
    )
    orig_dna = orig_sp + mature_skewed + "TAA"
    # Synthetic DP output with different codons in SP and skewed mature codons
    dp_sp = "ATGGCAAAGACTAACCTCTTCCTCTTCCTCATCTTCTCACTCCTACTCTCACTCTCATCAGCA"[:60]
    dp_cds = dp_sp + mature_skewed + "TAA"

    refined = apply_default_sop_refinement(
        dp_cds,
        original_dna=orig_dna,
        preserve_signal_peptide=True,
        gc_min=35.0,
        gc_max=55.0,
        forbidden_motifs=("GGTCTC", "GAAGAC", "CGTCTC"),
    )

    assert len(refined) == len(dp_cds)
    assert str(Seq(refined).translate()) == str(Seq(orig_dna).translate())
    assert refined[:60] == orig_sp
    assert refined.endswith("TAA")

    mature_codons = [refined[i : i + 3] for i in range(60, len(refined) - 3, 3)]
    arg_counts = {c: mature_codons.count(c) for c in ARG_CODONS}
    ser_counts = {c: mature_codons.count(c) for c in SER_CODONS}
    assert abs(arg_counts["AGA"] - arg_counts["AGG"]) <= 2 and arg_counts["AGG"] >= 2
    assert all(cnt >= 1 for cnt in ser_counts.values())

    # No adjacent identical Arg or Ser codons in mature region
    for i in range(len(mature_codons) - 1):
        if mature_codons[i] in ARG_CODONS or mature_codons[i] in SER_CODONS:
            assert mature_codons[i] != mature_codons[i + 1]


def test_protected_pause_positions_are_not_modified() -> None:
    # 25-codon sequence where index 21 (0-based) is protected
    cds = ("GCT" * 20) + "AGA" + "AGA" + "TCT" + "TCT" + "GCT"
    refined = apply_default_sop_refinement(
        cds,
        original_dna=None,
        preserve_signal_peptide=False,
        gc_min=30.0,
        gc_max=65.0,
        protected_codon_indices={20},
    )
    assert refined[20 * 3 : 21 * 3] == "AGA"
    assert str(Seq(refined).translate()) == str(Seq(cds).translate())


def test_mature_cds_eliminates_polya_are_and_splice_donor_motifs() -> None:
    # Build a valid mature CDS containing:
    # - ATTTA (AU-rich destabilizing): ATT + TAC -> ATTTAC
    # - AATAAA (premature PolyA): AAT + AAA
    # - AGGT AAGT (cryptic splice donor CAG|GTAAGT): CAG + GTA + AGT
    # - AAAAA (5-nt homopolymer): AAA + AAG
    sp = "ATG" + ("GCT" * 19)  # 60 nt
    mature = (
        "ATTTAC"      # ATTTA
        + "AATAAA"    # AATAAA
        + "ATTAAA"    # ATTAAA
        + "AATAAT"    # AATAAT
        + "CAGGTAAGT" # [AC]AGGT[AG]AGT
        + "AAAAAG"    # 5-nt A homopolymer
        + ("GCT" * 12)
    )
    cds = sp + mature + "TAA"
    refined = apply_default_sop_refinement(
        cds,
        original_dna=cds,
        preserve_signal_peptide=True,
        gc_min=35.0,
        gc_max=55.0,
        forbidden_motifs=("GGTCTC", "GAAGAC", "CGTCTC"),
    )
    mature_out = refined[60:-3]
    assert str(Seq(refined).translate()) == str(Seq(cds).translate())
    assert refined[:60] == sp
    for bad in ("ATTTA", "AATAAA", "ATTAAA", "AATAAT", "AAAAA"):
        assert bad not in mature_out
    import re
    assert not re.search(r"[AC]AGGT[AG]AGT", mature_out)
