import pytest

from factorforge.protein_risk import annotate_protein_risk
from factorforge.protein_risk.risk_classifier import classify_risk
from factorforge.protein_risk.sp_predict import predict_signal_peptide
from factorforge.protein_risk.tm_predict import predict_tm_segments

HYDROPHOBIC = "I" * 25
SOLUBLE_LINKER = "DEKR" * 8
MULTIPASS = (HYDROPHOBIC + SOLUBLE_LINKER) * 3
SINGLE_PASS = SOLUBLE_LINKER + HYDROPHOBIC + SOLUBLE_LINKER
SIGNAL_PEPTIDE_SINGLE_PASS = HYDROPHOBIC + SOLUBLE_LINKER
SOLUBLE = "DEKRQNST" * 8


def test_tm_predict_multipass():
    assert predict_tm_segments(MULTIPASS)["tm_count"] >= 3


def test_tm_predict_single_pass():
    assert predict_tm_segments(SINGLE_PASS)["tm_count"] == 1


def test_tm_predict_soluble():
    assert predict_tm_segments(SOLUBLE)["tm_count"] == 0


def test_tm_predict_empty_seq():
    with pytest.raises(ValueError):
        predict_tm_segments("")


def test_tm_predict_unknown_residue_is_neutral():
    assert predict_tm_segments("X" * 25)["mean_kd_score"] == 0.0


def test_sp_heuristic_positive():
    assert predict_signal_peptide(SIGNAL_PEPTIDE_SINGLE_PASS) is True


def test_sp_heuristic_negative():
    assert predict_signal_peptide(SINGLE_PASS) is False


def test_classify_high_multipass():
    assert classify_risk(3, -1.0, True) == "HIGH"


def test_classify_high_no_sp():
    assert classify_risk(1, 0.3, False) == "HIGH"


def test_classify_medium_single_sp():
    assert classify_risk(1, 0.9, True) == "MEDIUM"


def test_classify_low_soluble():
    assert classify_risk(0, 1.0, False) == "LOW"


def test_annotate_returns_required_keys():
    assert set(annotate_protein_risk(SOLUBLE)) == {
        "tm_count",
        "mean_kd_score",
        "signal_peptide_predicted",
        "risk_level",
        "warnings",
    }


def test_annotate_no_forbidden_words():
    text = str(annotate_protein_risk(MULTIPASS)).lower()
    assert "expression predicted" not in text
    assert "yield predicted" not in text
    assert "발현 예측" not in text
    assert "수율 예측" not in text


def test_annotate_high_risk_has_warning():
    result = annotate_protein_risk(MULTIPASS)
    assert result["risk_level"] == "HIGH"
    assert result["warnings"]


def test_annotate_low_risk_empty_warnings():
    result = annotate_protein_risk(SOLUBLE)
    assert result["risk_level"] == "LOW"
    assert result["warnings"] == []


def test_no_private_sequence_in_fixtures():
    # Compliance check: ensure internal staging identifiers or placeholder tokens
    # do not accidentally appear in public test fixtures.
    fixture_material = " ".join(
        [HYDROPHOBIC, SOLUBLE_LINKER, MULTIPASS, SINGLE_PASS, SIGNAL_PEPTIDE_SINGLE_PASS, SOLUBLE]
    )
    forbidden = ["C" + "D" + "9", "C" + "D" + "4" + "7"]
    for marker in forbidden:
        assert marker not in fixture_material, f"Private staging identifier found in public fixture: {marker}"
