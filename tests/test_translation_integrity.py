from factorforge.analysis.metrics import translate_dna
from factorforge.utils.sequence_validator import validate_cds_output


def test_translate_known():
    # ATG=M, AAA=K, ACC=T
    assert translate_dna("ATGAAAACC").rstrip("*") == "MKT"


def test_validate_passes_for_matching_cds():
    res = validate_cds_output("MKT", "ATGAAAACC")
    assert res["passed"] is True
    assert res["aa_identity"] == 1.0


def test_validate_fails_on_internal_stop():
    # TAA internal stop -> M*  (mismatch with MK)
    res = validate_cds_output("MK", "ATGTAAAAA")
    assert res["passed"] is False
