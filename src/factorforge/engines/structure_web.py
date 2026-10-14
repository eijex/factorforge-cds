"""Bounded, stateless structure intake for direct experimental web testing."""

from io import StringIO
import math
import warnings
from Bio.PDB import PDBParser, MMCIFParser
from Bio.SeqUtils import seq1
from factorforge.validation.tradeoff_suite import compute_tradeoff_comparison

MAX_STRUCTURE_BYTES = 250_000
MAX_STRUCTURE_PROTEIN_AA = 300


def run_structure_web(protein, context, weights, gc_min, gc_max, forbidden_motifs):
    if not isinstance(context, dict):
        raise ValueError("3DP requires a structure_context object")
    if not 3 <= len(protein) <= MAX_STRUCTURE_PROTEIN_AA:
        raise ValueError(
            "Direct 3DP testing supports 3-300 amino acids; use local reviewed workflows for larger inputs"
        )
    if context.get("positions_reviewed") is not True:
        raise ValueError("Confirm that the structure positions have been reviewed")
    text = context.get("content")
    if not isinstance(text, str) or not text or len(text.encode("utf-8")) > MAX_STRUCTURE_BYTES:
        raise ValueError("Upload a PDB/mmCIF structure of at most 250 KB")
    fmt = context.get("format")
    if fmt not in {"pdb", "mmcif"}:
        raise ValueError("Structure format must be pdb or mmcif")
    chain_id = context.get("chain")
    model_number = context.get("model", 1)
    if not isinstance(chain_id, str) or not 1 <= len(chain_id) <= 8:
        raise ValueError("An explicit structure chain is required")
    if type(model_number) is not int or model_number < 1:
        raise ValueError("Model must be a positive integer (one-based file model order)")
    positions = context.get("positions")
    if (
        not isinstance(positions, list)
        or not positions
        or any(type(p) is not int or not 1 <= p <= len(protein) for p in positions)
    ):
        raise ValueError(
            "Reviewed positions must be one-based protein positions within the input sequence"
        )
    if len(positions) != len(set(positions)):
        raise ValueError("Reviewed positions must be unique")
    target = context.get("target_adaptiveness", 0.40)
    penalty = context.get("penalty_weight", 2.0)
    if (
        isinstance(target, bool)
        or not isinstance(target, (int, float))
        or not math.isfinite(target)
        or not 0 < target <= 1
    ):
        raise ValueError("Target adaptiveness must be finite and in (0, 1]")
    if (
        isinstance(penalty, bool)
        or not isinstance(penalty, (int, float))
        or not math.isfinite(penalty)
        or not 0 <= penalty <= 100
    ):
        raise ValueError("Penalty weight must be finite and in [0, 100]")
    if not all(math.isfinite(x) for x in (gc_min, gc_max)) or not 0 <= gc_min <= gc_max <= 1:
        raise ValueError("3DP GC constraints must form a finite percentage range")
    try:
        parser = (
            PDBParser(PERMISSIVE=False, QUIET=True) if fmt == "pdb" else MMCIFParser(QUIET=True)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure("web-input", StringIO(text))
        models = list(structure.get_models())
        if model_number > len(models):
            raise ValueError()
        chain = models[model_number - 1][chain_id]
        residues = [r for r in chain.get_residues() if r.id[0] == " "]
        if any(r.is_disordered() == 2 for r in residues):
            raise ValueError()
        observed = "".join(seq1(r.resname, undef_code="X") for r in residues)
    except Exception:
        raise ValueError(
            "Cannot read the selected structure chain/model; check the uploaded file"
        ) from None
    if observed != protein:
        raise ValueError(
            "Selected structure chain must exactly match the full protein sequence; partial or homologous structures require reviewed mapping"
        )
    comparison = compute_tradeoff_comparison(
        protein,
        positions,
        weights,
        gc_min,
        gc_max,
        forbidden_motifs,
        lambda_weight=penalty,
        target_adaptiveness=target,
    )
    if not comparison.hard_constraints_pass:
        raise ValueError(
            "A generated candidate failed independent hard checks; no candidate was returned"
        )
    metadata = {
        "format": fmt,
        "chain": chain_id,
        "model": model_number,
        "mapped_residues": len(residues),
        "reviewed_positions": sorted(positions),
        "mapping_mode": "exact_full_chain",
        "target_adaptiveness": target,
        "penalty_weight": penalty,
    }
    return comparison, metadata
