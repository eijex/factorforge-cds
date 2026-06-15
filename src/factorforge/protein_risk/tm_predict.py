"""Deterministic transmembrane-segment heuristic."""

from factorforge.protein_risk.kd_scale import KD_SCALE, TM_THRESHOLD, TM_WINDOW


def predict_tm_segments(protein_seq: str) -> dict:
    """Return merged hydrophobic windows and whole-protein mean KD score."""
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    scores = [KD_SCALE.get(residue, 0.0) for residue in protein_seq.upper()]
    candidate_segments: list[tuple[int, int]] = []
    for start in range(max(0, len(scores) - TM_WINDOW + 1)):
        window_mean = sum(scores[start : start + TM_WINDOW]) / TM_WINDOW
        if window_mean >= TM_THRESHOLD:
            candidate_segments.append((start, start + TM_WINDOW))

    merged_segments: list[tuple[int, int]] = []
    for start, end in candidate_segments:
        if merged_segments and start <= merged_segments[-1][1]:
            previous_start, previous_end = merged_segments[-1]
            merged_segments[-1] = (previous_start, max(previous_end, end))
        else:
            merged_segments.append((start, end))

    return {
        "tm_count": len(merged_segments),
        "mean_kd_score": sum(scores) / len(scores),
        "segments": merged_segments,
    }
