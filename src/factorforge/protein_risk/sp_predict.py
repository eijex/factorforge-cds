"""Deterministic N-terminal signal-peptide heuristic."""

from factorforge.protein_risk.kd_scale import SP_HRUN_MIN, SP_SCAN_WINDOW

HYDROPHOBIC_RESIDUES = frozenset("ILVAMFW")


def predict_signal_peptide(protein_seq: str) -> bool:
    """Return whether an N-terminal hydrophobic run meets the heuristic."""
    run_length = 0
    for residue in protein_seq.upper()[:SP_SCAN_WINDOW]:
        if residue in HYDROPHOBIC_RESIDUES:
            run_length += 1
            if run_length >= SP_HRUN_MIN:
                return True
        else:
            run_length = 0
    return False
