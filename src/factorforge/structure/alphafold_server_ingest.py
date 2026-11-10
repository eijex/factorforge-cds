"""AlphaFold Server (alphafoldserver.com) structure ingest and confidence analyzer.

Parses AlphaFold 3 mmCIF (*_model_0.cif) output and summary confidence JSON,
extracts CA-atom pLDDT scores, classifies structural flexibility, and calculates
cotranslational ribosome exit-tunnel (+35 aa) pause recommendations.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Computational structural confidence · In-silico pause scheduling hypothesis
"""

from __future__ import annotations

import json
import math
from io import StringIO
from pathlib import Path
from typing import Any

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

DEFAULT_EXIT_TUNNEL_OFFSET = 35
PAUSE_MIN_REGION_LENGTH = 3
PAUSE_PLDDT_THRESHOLD = 70.0


class AlphaFoldIngestError(ValueError):
    """Raised when parsing or structural validation of AlphaFold Server output fails."""


def _read_local_text(source: str | Path) -> str:
    if isinstance(source, Path):
        return source.read_text(encoding="utf-8-sig")
    if len(source) <= 240 and "\n" not in source and not source.lstrip().startswith(("{", "[", "data_")):
        path = Path(source)
        if path.is_file():
            return path.read_text(encoding="utf-8-sig")
    return source.removeprefix("\ufeff")


def _pause_candidates(residues: list[dict[str, Any]], tunnel_offset: int) -> list[dict[str, Any]]:
    observed = {r["seq_id"] for r in residues}
    candidates = []
    run = []
    for residue in [*residues, None]:
        low_confidence = residue is not None and residue["plddt"] < PAUSE_PLDDT_THRESHOLD
        if run and (not low_confidence or residue["seq_id"] != run[-1]["seq_id"] + 1):
            if len(run) >= PAUSE_MIN_REGION_LENGTH:
                end = run[-1]["seq_id"]
                center = end + tunnel_offset
                candidates.append({
                    "start_residue": run[0]["seq_id"],
                    "end_residue": end,
                    "length": len(run),
                    "mean_plddt": round(sum(r["plddt"] for r in run) / len(run), 2),
                    "tunnel_pause_position": center,
                    "within_structure_chain": center in observed,
                    "within_precursor": None,
                })
            run = []
        if low_confidence:
            run.append(residue)
    return candidates


pause_candidates = _pause_candidates


def parse_af3_cif(

    cif_path_or_text: str | Path,
    target_chain: str = "A",
    tunnel_offset: int = DEFAULT_EXIT_TUNNEL_OFFSET,
) -> dict[str, Any]:
    """Parse AlphaFold 3 mmCIF structure file and extract CA-atom pLDDT confidence.

    Args:
        cif_path_or_text: File path or raw text string of the mmCIF file.
        target_chain: Chain identifier (e.g. 'A', 'B', or 'ALL').
        tunnel_offset: Ribosome exit-tunnel latency offset in amino acids (default: 35).

    Returns:
        Structured dictionary containing sequence, CA-pLDDT scores, distribution,
        and cotranslational pause recommendations.
    """
    if type(tunnel_offset) is not int or tunnel_offset < 0:
        raise AlphaFoldIngestError("The structure offset must be a nonnegative integer.")
    text = _read_local_text(cif_path_or_text)
    try:
        columns = MMCIF2Dict(StringIO(text))
    except ValueError as exc:
        raise AlphaFoldIngestError("Invalid mmCIF structure.") from exc

    atom = columns.get("_atom_site.label_atom_id")
    if not atom:
        raise AlphaFoldIngestError("No _atom_site loop found in mmCIF structure.")
    n_rows = len(atom)
    chains = columns.get("_atom_site.auth_asym_id", columns.get("_atom_site.label_asym_id"))
    positions = columns.get("_atom_site.label_seq_id")
    components = columns.get("_atom_site.label_comp_id")
    confidence = columns.get("_atom_site.B_iso_or_equiv")
    models = columns.get("_atom_site.pdbx_PDB_model_num", ["1"] * n_rows)
    groups = columns.get("_atom_site.group_PDB", ["ATOM"] * n_rows)
    if any(column is None or len(column) != n_rows for column in (
        chains, positions, components, confidence, models, groups,
    )):
        raise AlphaFoldIngestError("Missing or inconsistent mmCIF atom_site columns.")

    chains_data: dict[str, list[dict[str, Any]]] = {}
    seen: set[tuple[str, int]] = set()
    selected_models: set[str] = set()
    for name, chain, pos, comp, score, model, group in zip(
        atom, chains, positions, components, confidence, models, groups, strict=True,
    ):
        if name != "CA" or group != "ATOM" or (target_chain.upper() != "ALL" and chain != target_chain):
            continue
        try:
            seq_id = int(pos)
            plddt = float(score)
        except ValueError as exc:
            raise AlphaFoldIngestError("Invalid polymer position or CA confidence.") from exc
        if seq_id < 1 or not math.isfinite(plddt) or not 0 <= plddt <= 100:
            raise AlphaFoldIngestError("Polymer positions and CA confidence are out of range.")
        if comp not in THREE_TO_ONE:
            raise AlphaFoldIngestError("Unsupported polymer residue in selected chain.")
        selected_models.add(model)
        if (chain, seq_id) in seen or len(selected_models) != 1:
            raise AlphaFoldIngestError("Ambiguous CA mapping: duplicate residue or multiple models.")
        seen.add((chain, seq_id))
        chains_data.setdefault(chain, []).append({
            "seq_id": seq_id, "aa": THREE_TO_ONE[comp], "aa_3": comp, "plddt": plddt,
        })

    if not chains_data:
        raise AlphaFoldIngestError(f"No CA atoms found for target chain '{target_chain}'.")

    # Analyze chains
    analyzed_chains: dict[str, Any] = {}

    for ch_name, residues in chains_data.items():
        residues.sort(key=lambda r: r["seq_id"])
        seq = "".join(r["aa"] for r in residues)
        plddts = [r["plddt"] for r in residues]
        length = len(residues)
        mean_plddt = round(sum(plddts) / length, 2) if length > 0 else 0.0

        v_high = sum(1 for p in plddts if p >= 90.0)
        confident = sum(1 for p in plddts if 70.0 <= p < 90.0)
        low = sum(1 for p in plddts if 50.0 <= p < 70.0)
        v_low = sum(1 for p in plddts if p < 50.0)

        candidates = _pause_candidates(residues, tunnel_offset)

        analyzed_chains[ch_name] = {
            "chain_id": ch_name,
            "length": length,
            "sequence": seq,
            "mean_plddt": mean_plddt,
            "plddt_distribution": {
                "very_high_ge_90": {"count": v_high, "percent": round(v_high / length * 100, 1)},
                "confident_70_89": {"count": confident, "percent": round(confident / length * 100, 1)},
                "low_50_69": {"count": low, "percent": round(low / length * 100, 1)},
                "very_low_lt_50": {"count": v_low, "percent": round(v_low / length * 100, 1)},
            },
            "pause_candidates": candidates,
            "residues": residues,
        }

    # If single chain requested and found, surface as primary
    primary_chain = target_chain if target_chain.upper() != "ALL" and target_chain in analyzed_chains else next(iter(analyzed_chains.keys()))

    return {
        "status": "AVAILABLE",
        "primary_chain": primary_chain,
        "chains": analyzed_chains,
        "chain_count": len(analyzed_chains),
        "tunnel_offset_applied": tunnel_offset,
        "coordinate_system": "mmCIF_label_seq_id_1based",
        "precursor_mapping": "NOT_VERIFIED",
        "epistemic_status": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
    }


def parse_af3_summary_confidences(json_path_or_text: str | Path) -> dict[str, Any]:
    """Parse AlphaFold 3 summary_confidences_0.json file.

    Extracts global ranking_score, ptm, iptm, and disorder fraction.
    """
    raw_text = _read_local_text(json_path_or_text)
    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError as exc:
        raise AlphaFoldIngestError("Invalid summary confidences JSON.") from exc
    if not isinstance(data, dict):
        raise AlphaFoldIngestError("Summary confidences must be a JSON object.")

    metrics = {}
    for key in ("ranking_score", "ptm", "iptm", "fraction_disordered"):
        value = data.get(key)
        if value is not None and (
            type(value) not in (int, float) or not math.isfinite(value)
            or (key != "ranking_score" and not 0 <= value <= 1)
        ):
            raise AlphaFoldIngestError(f"Invalid summary confidence field: {key}.")
        metrics[key] = value
    clash = data.get("has_clash")
    if clash is not None and (type(clash) not in (int, float, bool) or clash not in (0, 1)):
        raise AlphaFoldIngestError("Invalid summary confidence field: has_clash.")
    return {**metrics, "has_clash": None if clash is None else bool(clash)}


def ingest_alphafold_server_result(
    cif_path: Path | str,
    chain: str = "A",
    summary_json_path: Path | str | None = None,
    output_path: Path | str | None = None,
    tunnel_offset: int = DEFAULT_EXIT_TUNNEL_OFFSET,
) -> dict[str, Any]:
    """Ingest AlphaFold Server prediction results and generate structured sidecar."""
    structure_data = parse_af3_cif(Path(cif_path), target_chain=chain, tunnel_offset=tunnel_offset)

    summary_metrics: dict[str, Any] | None = None
    if summary_json_path:
        summary_metrics = parse_af3_summary_confidences(Path(summary_json_path))

    result = {
        **structure_data,
        "summary_confidences": summary_metrics,
        "summary_binding": "UNVERIFIED" if summary_metrics is not None else "NOT_PROVIDED",
    }

    if output_path:
        out_target = Path(output_path)
        out_target.parent.mkdir(parents=True, exist_ok=True)
        with out_target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")

    return result


def format_ingest_report(result: dict[str, Any]) -> str:
    """Format the ingested AlphaFold Server result into clear human-readable terminal text."""
    lines: list[str] = [
        "=== AlphaFold Server (AF3) Ingest Report ===",
        "Coordinate system: model polymer numbering. Precursor mapping not verified.",
        f"Chains Analyzed: {result.get('chain_count')}",
    ]

    summary = result.get("summary_confidences")
    if summary:
        lines.append(
            f"Unbound summary · Ranking Score: {summary.get('ranking_score')} | "
            f"pTM: {summary.get('ptm')} | "
            f"ipTM: {summary.get('iptm')}"
        )

    chains = result.get("chains", {})
    for ch_id, ch in chains.items():
        lines.append("")
        lines.append(f"Chain '{ch_id}' (Length: {ch['length']} aa, Mean pLDDT: {ch['mean_plddt']}):")
        dist = ch["plddt_distribution"]
        lines.append(
            f"  pLDDT Bands: Very High(>=90): {dist['very_high_ge_90']['percent']}% | "
            f"Confident(70-89): {dist['confident_70_89']['percent']}% | "
            f"Low(50-69): {dist['low_50_69']['percent']}% | "
            f"Very Low(<50): {dist['very_low_lt_50']['percent']}%"
        )

        candidates = ch.get("pause_candidates", [])
        if candidates:
            lines.append(f"  Offset Hypotheses (+{result['tunnel_offset_applied']} aa exit-tunnel offset):")
            for c in candidates:
                status_str = "observed model position" if c["within_structure_chain"] else "position not observed"
                lines.append(
                    f"    - Residues {c['start_residue']}-{c['end_residue']} (pLDDT: {c['mean_plddt']}) "
                    f"-> Model Offset Position: {c['tunnel_pause_position']} ({status_str})"
                )
        else:
            lines.append("  Offset Hypotheses: None in observed residues")

    lines.append("")
    lines.append(f"Epistemic: {result['epistemic_status']} · Biological validation: {result['biological_validation']}")
    return "\n".join(lines)


def extract_pause_positions_from_sidecar(
    sidecar_path_or_data: str | Path | dict[str, Any],
    protein_sequence: str,
    allow_mock: bool = False,
) -> list[int]:
    """Bind a single, complete model chain to the input and derive eligible offsets."""
    if isinstance(sidecar_path_or_data, (str, Path)):
        try:
            data = json.loads(Path(sidecar_path_or_data).read_text(encoding="utf-8-sig"))
        except FileNotFoundError as exc:
            raise AlphaFoldIngestError("Structure sidecar file not found.") from exc
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise AlphaFoldIngestError("Failed to read structure sidecar JSON.") from exc
    else:
        data = sidecar_path_or_data

    if not isinstance(data, dict) or data.get("status") != "AVAILABLE":
        raise AlphaFoldIngestError("Sidecar status is not AVAILABLE or structure is invalid.")
    if (data.get("is_mock") is True or data.get("source_type") == "SYNTHETIC_FIXTURE") and not allow_mock:
        raise AlphaFoldIngestError(
            "Sidecar was generated from synthetic/mock data. "
            "Mock sidecars are blocked from production/research runs. "
            "Pass --allow-mock to proceed with synthetic testing."
        )
    if data.get("coordinate_system") not in ("mmCIF_label_seq_id_1based", "PDB_res_seq_1based"):
        raise AlphaFoldIngestError("Unsupported coordinate system. Regenerate the sidecar.")
    chains = data.get("chains")
    chain_id = data.get("primary_chain")
    if not isinstance(chains, dict) or len(chains) != 1 or not isinstance(chain_id, str):
        raise AlphaFoldIngestError("Sidecar must contain exactly one selected chain. Use ingest-af3 --chain.")
    chain = chains.get(chain_id)
    protein = "".join(protein_sequence.split()).upper().removesuffix("*")
    if not protein or set(protein) - set(THREE_TO_ONE.values()):
        raise AlphaFoldIngestError("The optimization input must be a standard amino acid sequence.")
    if not isinstance(chain, dict) or chain.get("chain_id") != chain_id or chain.get("sequence") != protein:
        raise AlphaFoldIngestError("Sidecar chain sequence does not match the optimization input.")
    residues = chain.get("residues")
    if not isinstance(residues, list) or len(residues) != len(protein) or chain.get("length") != len(protein):
        raise AlphaFoldIngestError("A complete one-to-one model residue mapping is required.")
    for index, (residue, aa) in enumerate(zip(residues, protein, strict=True), start=1):
        if (
            not isinstance(residue, dict)
            or type(residue.get("seq_id")) is not int
            or residue["seq_id"] != index
            or residue.get("aa") != aa
        ):
            raise AlphaFoldIngestError("Model residues must match the input at consecutive one-based positions.")
        score = residue.get("plddt")
        if type(score) not in (int, float) or not 0 <= score <= 100:
            raise AlphaFoldIngestError("Invalid CA confidence in sidecar residue mapping.")
    offset = data.get("tunnel_offset_applied")
    if type(offset) is not int or offset < 0:
        raise AlphaFoldIngestError("Sidecar offset must be a nonnegative integer.")

    candidates = _pause_candidates(residues, offset)
    if chain.get("pause_candidates") != candidates:
        raise AlphaFoldIngestError("Sidecar hypotheses disagree with the residue data. Regenerate with ingest-af3.")
    positions = sorted({c["tunnel_pause_position"] for c in candidates if c["within_structure_chain"]})
    if not positions:
        raise AlphaFoldIngestError("No in-bounds model offset positions. Use plain DP for a baseline.")
    return positions
