"""Compose protein structural-risk annotations."""

from factorforge.protein_risk.risk_classifier import classify_risk
from factorforge.protein_risk.sp_predict import predict_signal_peptide
from factorforge.protein_risk.tm_predict import predict_tm_segments


def annotate_protein_risk(protein_seq: str) -> dict:
    """Return deterministic structural-risk indicators without sequence material."""
    tm_prediction = predict_tm_segments(protein_seq)
    sp_predicted = predict_signal_peptide(protein_seq)
    risk_level = classify_risk(
        tm_prediction["tm_count"],
        tm_prediction["mean_kd_score"],
        sp_predicted,
    )

    warnings = []
    if risk_level != "LOW":
        warnings.append(
            "Protein structural risk indicator: "
            f"{risk_level} (tm_count={tm_prediction['tm_count']}). "
            "CDS-level checks are independent of protein structural risk."
        )

    return {
        "tm_count": tm_prediction["tm_count"],
        "mean_kd_score": tm_prediction["mean_kd_score"],
        "signal_peptide_predicted": sp_predicted,
        "risk_level": risk_level,
        "warnings": warnings,
    }
