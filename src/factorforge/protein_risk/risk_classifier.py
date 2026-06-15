"""Protein structural-risk classification rules."""

from factorforge.protein_risk.kd_scale import MEAN_KD_RISK_THRESHOLD


def classify_risk(tm_count: int, mean_kd: float, sp_predicted: bool) -> str:
    """Classify deterministic structural-risk indicators."""
    if tm_count >= 3 or (
        tm_count >= 1 and mean_kd >= MEAN_KD_RISK_THRESHOLD and not sp_predicted
    ):
        return "HIGH"
    if tm_count >= 1:
        return "MEDIUM"
    return "LOW"
