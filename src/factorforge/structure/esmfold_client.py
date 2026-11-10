"""ESMFold structure prediction client and sidecar generator.

Queries a configured ESMFold service, extracts CA confidence, and generates
position-offset hypotheses.

EPISTEMIC NOTICE:
UNVERIFIED_PUBLIC_CORPUS · Biological validation: NOT_AVAILABLE
Computational neural structure prediction · In-silico pause scheduling hypothesis
Low pLDDT (<70) reflects computational model uncertainty or predicted structural
disorder, not confirmed biological hinges or in-vivo ribosomal dwell kinetics.
Explicit wet-lab validation is required before drawing biological conclusions.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import requests

from factorforge.structure.alphafold_server_ingest import (
    DEFAULT_EXIT_TUNNEL_OFFSET,
    PAUSE_PLDDT_THRESHOLD,
    THREE_TO_ONE,
    pause_candidates,
)

ESM_CLIENT_VERSION = "1.0.0"
ESM_ATLAS_DEFAULT_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}


class ESMFoldClientError(ValueError):
    """Raised when ESMFold API communication or PDB parsing fails."""


def _clean_sequence(sequence: str) -> str:
    """Normalize amino acid sequence and validate standard alphabet."""
    cleaned = "".join(sequence.split()).upper().removesuffix("*")
    if not cleaned:
        raise ESMFoldClientError("Protein sequence cannot be empty.")
    invalid = set(cleaned) - set(THREE_TO_ONE.values())
    if invalid:
        raise ESMFoldClientError(f"Sequence contains non-standard amino acid codes: {sorted(invalid)}")
    return cleaned


def _generate_mock_pdb(sequence: str, chain: str = "A") -> str:
    """Generate deterministic mock PDB for testing and offline environments.

    Assigns high confidence (85.0) to structured regions and low confidence (45.0)
    to a flexible loop around residue 1..3 (or 10..15 depending on length) to ensure
    pause candidates can be exercised deterministically.
    """
    clean_seq = _clean_sequence(sequence)
    lines = ["HEADER    MOCK ESMFOLD PREDICTION", f"COMPND    CHAIN {chain}"]
    atom_id = 1

    # Insert flexible loop at residue 1..3 if length >= 10, else middle
    for i, aa in enumerate(clean_seq, start=1):
        res3 = ONE_TO_THREE.get(aa, "ALA")
        # Flexible region for testing: first 3 residues or residue 10-13
        plddt = 45.0 if (1 <= i <= 3 or 10 <= i <= 13) else 88.5

        # Mock synthetic alpha helix / extended coordinates
        x = round(math.cos(i * 0.5) * 5.0, 3)
        y = round(math.sin(i * 0.5) * 5.0, 3)
        z = round(i * 1.5, 3)

        # Standard PDB ATOM record format (80 columns)
        line = (
            f"ATOM  {atom_id:5d}  CA  {res3:>3s} {chain}{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{plddt:6.2f}           C"
        )
        lines.append(line)
        atom_id += 1

    lines.append("END")
    return "\n".join(lines) + "\n"


def predict_structure_esm(
    sequence: str,
    provider: str = "esmatlas",
    api_key: str | None = None,
    timeout: float = 30.0,
) -> str:
    """Query ESMFold API to obtain 3D structure in PDB format.

    Args:
        sequence: Amino acid sequence string.
        provider: 'esmatlas', 'nvidia', or 'mock'.
        api_key: Optional API key for NVIDIA BioNeMo (defaults to env vars).
        timeout: Network timeout in seconds.

    Returns:
        PDB text format string.
    """
    clean_seq = _clean_sequence(sequence)

    if provider == "mock":
        return _generate_mock_pdb(clean_seq)

    if not math.isfinite(timeout) or timeout <= 0:
        raise ESMFoldClientError("Timeout must be a finite positive number.")

    headers = {"User-Agent": f"FactorForge/{ESM_CLIENT_VERSION} (esmfold-client)"}
    if provider == "esmatlas":
        url = os.environ.get("ESM_ATLAS_URL", ESM_ATLAS_DEFAULT_URL)
        headers["Content-Type"] = "text/plain"
        body = clean_seq
    elif provider == "nvidia":
        url = os.environ.get("NVIDIA_ESMFOLD_URL", "")
        if not url:
            raise ESMFoldClientError("Set NVIDIA_ESMFOLD_URL to a configured ESMFold service.")
        token = api_key or os.environ.get("NVCF_RUN_KEY") or os.environ.get("NVIDIA_API_KEY")
        if not token:
            raise ESMFoldClientError("Set NVCF_RUN_KEY or NVIDIA_API_KEY, or pass --api-key.")
        headers.update(Authorization=f"Bearer {token}")
        headers["Content-Type"] = "application/json"
        body = json.dumps({"sequence": clean_seq})
    else:
        raise ESMFoldClientError(f"Unknown provider '{provider}'. Choose esmatlas, nvidia, or mock.")

    try:
        endpoint = urlsplit(url)
        if endpoint.scheme != "https" or not endpoint.hostname or endpoint.username is not None:
            raise ESMFoldClientError("The remote endpoint must use HTTPS without URL credentials.")
        response = requests.post(
            url, data=body, headers=headers, timeout=timeout, allow_redirects=False
        )
        if response.status_code != 200:
            raise ESMFoldClientError(f"{provider} returned HTTP {response.status_code}.")
        if provider == "esmatlas":
            pdb = response.text
        else:
            payload = response.json()
            pdb = payload.get("pdb") if isinstance(payload, dict) else None
            if isinstance(payload, dict) and isinstance(payload.get("pdbs"), list) and payload["pdbs"]:
                pdb = payload["pdbs"][0]
    except requests.RequestException as exc:
        raise ESMFoldClientError(f"{provider} request failed ({type(exc).__name__}).") from exc
    if not isinstance(pdb, str) or not pdb.strip():
        raise ESMFoldClientError(f"{provider} did not return PDB text.")
    return pdb


def parse_esm_pdb(
    pdb_text: str,
    target_chain: str = "A",
    tunnel_offset: int = DEFAULT_EXIT_TUNNEL_OFFSET,
    plddt_threshold: float = PAUSE_PLDDT_THRESHOLD,
    expected_sequence: str | None = None,
) -> dict[str, Any]:
    """Parse ESMFold PDB text and extract CA-pLDDT confidence metrics.

    In ESMFold PDBs, per-residue pLDDT scores (0 to 100) are stored in the
    temperature factor (B-factor) column.
    """
    if type(tunnel_offset) is not int or tunnel_offset < 0:
        raise ESMFoldClientError("Tunnel offset must be a nonnegative integer.")
    if plddt_threshold != PAUSE_PLDDT_THRESHOLD:
        raise ESMFoldClientError("Sidecar requires the shared CA confidence threshold.")
    if not target_chain or target_chain == "ALL":
        raise ESMFoldClientError("Select one PDB chain.")
    if not pdb_text or not pdb_text.strip():
        raise ESMFoldClientError("PDB content is empty.")

    residues: list[dict[str, Any]] = []
    seen_seq_ids: set[int] = set()

    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 66:
            continue

        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        chain_id = line[21].strip() or "A"
        if chain_id != target_chain:
            continue

        try:
            seq_id = int(line[22:26].strip())
            res_name = line[17:20].strip()
            plddt = float(line[60:66].strip())
        except (ValueError, IndexError) as exc:
            raise ESMFoldClientError("Malformed PDB CA atom record.") from exc

        if math.isnan(plddt) or math.isinf(plddt) or not 0 <= plddt <= 100:
            raise ESMFoldClientError("CA confidence must be finite and between 0 and 100.")
        if line[26].strip() or seq_id in seen_seq_ids:
            raise ESMFoldClientError("Ambiguous CA residue numbering in the selected chain.")
        seen_seq_ids.add(seq_id)

        aa_1 = THREE_TO_ONE.get(res_name.upper())
        if aa_1 is None:
            raise ESMFoldClientError("Nonstandard CA residue in the selected chain.")
        residues.append({
            "seq_id": seq_id,
            "aa": aa_1,
            "plddt": round(plddt, 2),
        })

    if not residues:
        raise ESMFoldClientError(f"No CA atoms found for target chain '{target_chain}' in PDB.")

    # Sort residues by seq_id
    residues.sort(key=lambda r: r["seq_id"])

    # Ensure 1-based contiguous indices
    for expected_idx, r in enumerate(residues, start=1):
        if r["seq_id"] != expected_idx:
            raise ESMFoldClientError(
                f"Residue indices must be 1-based contiguous. Expected {expected_idx}, got {r['seq_id']}."
            )

    sequence = "".join(r["aa"] for r in residues)
    if expected_sequence is not None:
        clean_exp = _clean_sequence(expected_sequence)
        if sequence != clean_exp:
            raise ESMFoldClientError(
                f"Parsed PDB sequence does not match expected sequence "
                f"(length {len(sequence)} vs {len(clean_exp)})."
            )
    plddt_values = [r["plddt"] for r in residues]
    mean_plddt = round(sum(plddt_values) / len(plddt_values), 2)
    min_plddt = round(min(plddt_values), 2)
    max_plddt = round(max(plddt_values), 2)

    # Compute pause candidates using shared cotranslational tunnel offset logic
    candidates = pause_candidates(residues, tunnel_offset)

    primary_chain = target_chain
    return {
        "status": "AVAILABLE",
        "source_type": "ESMFold",
        "tool_version": ESM_CLIENT_VERSION,
        "coordinate_system": "PDB_res_seq_1based",
        "primary_chain": primary_chain,
        "tunnel_offset_applied": tunnel_offset,
        "plddt_threshold": plddt_threshold,
        "epistemic_status": "UNVERIFIED_PUBLIC_CORPUS",
        "biological_validation": "NOT_AVAILABLE",
        "reviewed": False,
        "chains": {
            primary_chain: {
                "chain_id": primary_chain,
                "length": len(residues),
                "sequence": sequence,
                "residues": residues,
                "pause_candidates": candidates,
                "metrics": {
                    "mean_plddt": mean_plddt,
                    "min_plddt": min_plddt,
                    "max_plddt": max_plddt,
                },
            }
        },
    }


def predict_and_ingest_esm(
    sequence: str,
    provider: str = "esmatlas",
    api_key: str | None = None,
    tunnel_offset: int = DEFAULT_EXIT_TUNNEL_OFFSET,
    output_sidecar_path: Path | None = None,
    output_pdb_path: Path | None = None,
    timeout: float = 30.0,
) -> tuple[dict[str, Any], str]:
    """Execute end-to-end ESMFold prediction, parse PDB, and generate sidecar."""
    outputs = [p for p in (output_pdb_path, output_sidecar_path) if p is not None]
    if len({p.resolve() for p in outputs}) != len(outputs):
        raise ESMFoldClientError("Choose different output paths for the PDB and sidecar.")
    if any(p.exists() for p in outputs):
        raise ESMFoldClientError("Output already exists. Choose a different path.")

    protein = _clean_sequence(sequence)
    pdb_text = predict_structure_esm(
        sequence=protein,
        provider=provider,
        api_key=api_key,
        timeout=timeout,
    )

    result = parse_esm_pdb(
        pdb_text=pdb_text,
        target_chain="A",
        tunnel_offset=tunnel_offset,
        expected_sequence=protein,
    )
    result["provider"] = provider
    result["is_mock"] = (provider == "mock")
    if provider == "mock":
        result["source_type"] = "SYNTHETIC_FIXTURE"
        result["provenance_basis"] = "MOCK_SYNTHETIC_DATA"

    created: list[Path] = []
    try:
        for path, content in (
            (output_pdb_path, pdb_text),
            (output_sidecar_path, json.dumps(result, indent=2, allow_nan=False)),
        ):
            if path is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("x", encoding="utf-8", newline="") as stream:
                    created.append(path)
                    stream.write(content)
    except OSError:
        for path in created:
            path.unlink(missing_ok=True)
        raise

    return result, pdb_text


def format_esm_report(result: dict[str, Any]) -> str:
    """Format human-readable CLI summary report for ESMFold prediction."""
    primary = result.get("primary_chain", "A")
    chain = result.get("chains", {}).get(primary, {})
    metrics = chain.get("metrics", {})
    candidates = chain.get("pause_candidates", [])

    lines = [
        "=======================================================",
        f"ESMFold Structure Report (v{result.get('tool_version', '1.0.0')})",
        f"Provider: {result.get('provider', 'esmatlas')} · Chain: {primary}",
        "=======================================================",
        f"Protein Length: {chain.get('length', 0)} aa",
        f"Mean CA pLDDT: {metrics.get('mean_plddt', 0.0):.2f} (Min: {metrics.get('min_plddt', 0.0):.2f}, Max: {metrics.get('max_plddt', 0.0):.2f})",
        f"Ribosome Tunnel Offset Applied: +{result.get('tunnel_offset_applied', 35)} aa",
        "",
    ]

    if result.get("is_mock"):
        lines.append("[WARNING: MOCK / SYNTHETIC DATA - NOT SUITABLE FOR PRODUCTION OR BIOLOGICAL DESIGNS]")
        lines.append("To use with factorforge optimize, you must explicitly pass --allow-mock.")
        lines.append("")
    if candidates:
        lines.append(f"Cotranslational Pause Hypotheses ({len(candidates)} regions):")
        for c in candidates:
            status = "in-bounds model position" if c.get("within_structure_chain") else "out-of-bounds"
            lines.append(
                f"  - Residues {c['start_residue']}-{c['end_residue']} "
                f"(Mean pLDDT: {c['mean_plddt']:.1f}) -> Tunnel Pause Position: {c['tunnel_pause_position']} ({status})"
            )
    else:
        lines.append("Cotranslational Pause Hypotheses: None under the CA confidence rule.")

    lines.append("")
    lines.append(f"Epistemic: {result['epistemic_status']} · Biological validation: {result['biological_validation']}")
    lines.append("Human-in-the-loop: reviewed=False (explicit review required before synthesis)")
    return "\n".join(lines)
