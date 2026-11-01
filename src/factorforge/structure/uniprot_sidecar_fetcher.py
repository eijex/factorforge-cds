"""UniProt Automated Sidecar Fetcher & Feature Resolver for FactorForge.

Fetches and normalizes experimental and curated protein annotations from UniProtKB
REST API into the exact sidecar schema required by `factorforge.nbfold.annotations`
and 3DP/4DP structural pause conditioning.

Epistemic & Governance Constraints:
- Epistemic classification: UNVERIFIED_PUBLIC_CORPUS / CURATED_UNIPROT_RECORD.
- Preserves explicit evidence codes (e.g. ECO:0000269 experimental vs ECO:0000255 sequence model).
- Clear boundary: UniProt feature presence describes protein biology, but translational pause
  site recommendation is an unverified computational heuristic (Biological validation: NOT_AVAILABLE).
- Mandatory human-in-the-loop: all fetched profiles initialize with reviewed=False.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Any, Dict, List, Optional
import urllib.request
import urllib.error

from factorforge.structure.auto_annotator import CandidateRegion

UNIPROT_FETCHER_VERSION = "1.1.0"
UNIPROT_API_ENDPOINT = "https://rest.uniprot.org/uniprotkb/{accession}.json"
ACCESSION_REGEX = re.compile(r"^[A-Z0-9_-]{4,20}$", re.IGNORECASE)


def fetch_uniprot_raw_record(
    accession: str,
    timeout: float = 8.0,
    headers: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Fetch raw JSON record directly from UniProt REST API."""
    clean_acc = accession.strip().upper()
    if not clean_acc or not ACCESSION_REGEX.match(clean_acc):
        raise ValueError(f"Invalid UniProt accession format: '{accession}'")

    url = UNIPROT_API_ENDPOINT.format(accession=clean_acc)
    req_headers = {"User-Agent": "FactorForge-UniProtFetcher/1.1 (Eijex; scientific-research)"}
    if headers:
        req_headers.update(headers)

    req = urllib.request.Request(url, headers=req_headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        content = resp.read().decode("utf-8")
        return json.loads(content)


class UniProtSidecarFetcher:
    """Fetches, caches, validates, and formats UniProt records into valid FactorForge sidecar contracts."""

    def __init__(self, cache_dir: Optional[Path] = None, timeout: float = 6.0) -> None:
        self.cache_dir = cache_dir
        self.timeout = timeout
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def fetch_and_format_sidecar(
        self,
        accession: str,
        target_protein_sequence: Optional[str] = None,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Fetch UniProt entry, perform strict validation, and format into sidecar payload."""
        if not isinstance(accession, str):
            return {
                "status": "ERROR_INVALID_ACCESSION",
                "reason": "Accession must be a non-empty string.",
                "accession": str(accession),
                "features": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        clean_acc = accession.strip().upper()
        if not clean_acc or not ACCESSION_REGEX.match(clean_acc):
            return {
                "status": "ERROR_INVALID_ACCESSION",
                "reason": f"Invalid UniProt accession format: '{accession}'",
                "accession": clean_acc,
                "features": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        # Input sequence validation if provided
        if target_protein_sequence is not None:
            if not isinstance(target_protein_sequence, str) or not target_protein_sequence.strip():
                return {
                    "status": "ERROR_INVALID_SEQUENCE",
                    "reason": "target_protein_sequence must be a non-empty string.",
                    "accession": clean_acc,
                    "features": [],
                    "candidate_regions": [],
                    "recommended_positions": [],
                    "reviewed": False,
                }

        cache_file = self.cache_dir / f"{clean_acc}_uniprot.json" if self.cache_dir else None
        raw_data: Optional[Dict[str, Any]] = None

        # Cache check with strict accession and JSON validation
        if use_cache and cache_file and cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    loaded = json.load(f)
                cached_acc = str(loaded.get("primaryAccession", "")).strip().upper()
                if cached_acc == clean_acc:
                    raw_data = loaded
                else:
                    # Mismatched cache content: invalidate
                    raw_data = None
            except Exception:
                raw_data = None

        if raw_data is None:
            try:
                raw_data = fetch_uniprot_raw_record(clean_acc, timeout=self.timeout)
                if cache_file and raw_data:
                    with open(cache_file, "w", encoding="utf-8") as f:
                        json.dump(raw_data, f, indent=2)
            except Exception as exc:
                return {
                    "status": "UNAVAILABLE",
                    "reason": f"Failed to fetch UniProt entry '{clean_acc}': {exc}",
                    "accession": clean_acc,
                    "features": [],
                    "candidate_regions": [],
                    "recommended_positions": [],
                    "reviewed": False,
                }

        seq_info = raw_data.get("sequence", {})
        seq_val = seq_info.get("value", "") if isinstance(seq_info, dict) else ""
        seq_clean = "".join(seq_val.upper().split()).rstrip("*")

        if not seq_clean:
            return {
                "status": "ERROR_EMPTY_RECORD_SEQUENCE",
                "reason": f"UniProt entry '{clean_acc}' contains no valid sequence.",
                "accession": clean_acc,
                "features": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        # Verify against input protein sequence
        if target_protein_sequence is not None:
            clean_input = "".join(target_protein_sequence.upper().split()).rstrip("*")
            if clean_input != seq_clean:
                return {
                    "status": "ERROR_SEQUENCE_MISMATCH",
                    "reason": (
                        f"UniProt sequence for '{clean_acc}' (len {len(seq_clean)}) "
                        f"does not match input target sequence (len {len(clean_input)})."
                    ),
                    "accession": clean_acc,
                    "features": [],
                    "candidate_regions": [],
                    "recommended_positions": [],
                    "reviewed": False,
                }

        digest = hashlib.sha256(seq_clean.encode("utf-8")).hexdigest()

        # Audit version handling: preserve exact values or None (never default missing version to 1)
        raw_audit = raw_data.get("entryAudit")
        entry_audit: Optional[Dict[str, Any]] = None
        if isinstance(raw_audit, dict):
            entry_audit = {
                "entryVersion": raw_audit.get("entryVersion"),
                "sequenceVersion": raw_audit.get("sequenceVersion"),
            }

        raw_features = raw_data.get("features", [])
        valid_features: List[Dict[str, Any]] = []
        candidate_regions: List[CandidateRegion] = []
        recommended_positions: List[int] = []

        for feat in raw_features:
            if not isinstance(feat, dict):
                continue
            ftype = feat.get("type", "")
            loc = feat.get("location", {})
            start_val: Optional[int] = None
            end_val: Optional[int] = None

            if isinstance(loc, dict):
                s_obj = loc.get("start")
                e_obj = loc.get("end")
                if isinstance(s_obj, dict):
                    start_val = s_obj.get("value")
                elif isinstance(s_obj, int):
                    start_val = s_obj
                if isinstance(e_obj, dict):
                    end_val = e_obj.get("value")
                elif isinstance(e_obj, int):
                    end_val = e_obj

            # Strict coordinate validation: integers, 1-based, within sequence bounds, start <= end
            if (
                isinstance(start_val, int)
                and isinstance(end_val, int)
                and 1 <= start_val <= end_val <= len(seq_clean)
            ):
                valid_features.append(feat)
                desc = feat.get("description", "")

                # Extract evidence code if present
                ev_code = "UNKNOWN"
                evidences = feat.get("evidences", [])
                if isinstance(evidences, list) and evidences:
                    first_ev = evidences[0]
                    if isinstance(first_ev, dict):
                        ev_code = first_ev.get("evidenceCode", "UNKNOWN")

                # Heuristic pause proposals based on protein features
                # Note: feature existence does NOT validate pause benefit (Biological validation: NOT_AVAILABLE)
                if ftype == "Signal":
                    recommended_positions.append(end_val)
                    candidate_regions.append(
                        CandidateRegion(
                            start_1based=start_val,
                            end_1based=end_val,
                            length=end_val - start_val + 1,
                            region_type="uniprot_signal_peptide",
                            reason=f"UniProt Signal peptide cleavage at {end_val} ({desc or 'cleavage site'}; {ev_code})",
                            heuristic_basis="UNIPROT_FEATURE_HEURISTIC_PROPOSAL",
                        )
                    )
                elif ftype == "Transmembrane":
                    recommended_positions.extend([start_val, end_val])
                    candidate_regions.append(
                        CandidateRegion(
                            start_1based=start_val,
                            end_1based=end_val,
                            length=end_val - start_val + 1,
                            region_type="uniprot_transmembrane_domain",
                            reason=f"UniProt Transmembrane helix ({desc or f'{start_val}-{end_val}'}; {ev_code})",
                            heuristic_basis="UNIPROT_FEATURE_HEURISTIC_PROPOSAL",
                        )
                    )
                elif ftype == "Disulfide bond":
                    recommended_positions.append(end_val)
                    candidate_regions.append(
                        CandidateRegion(
                            start_1based=start_val,
                            end_1based=end_val,
                            length=end_val - start_val + 1,
                            region_type="uniprot_disulfide_bond",
                            reason=f"UniProt Intrachain disulfide bridge (Cys{start_val}-Cys{end_val}; {ev_code})",
                            heuristic_basis="UNIPROT_FEATURE_HEURISTIC_PROPOSAL",
                        )
                    )

        unique_positions = sorted(list(set(recommended_positions)))

        # Build sidecar compliant with both nbfold and FactorForge structure pipelines
        sidecar_data: Dict[str, Any] = {
            "primaryAccession": raw_data.get("primaryAccession", clean_acc),
            "accession": clean_acc,
            "entryAudit": entry_audit,
            "sequence": {"value": seq_clean, "length": len(seq_clean)},
            "protein_sha256": digest,
            "protein_length": len(seq_clean),
            "features": valid_features,
            "uniprot_features": valid_features,  # Explicit descriptive biological feature table
            "pause_recommendation_candidates": [asdict_cand(cr) for cr in candidate_regions],  # Explicit heuristic candidates
            "candidate_regions": [asdict_cand(cr) for cr in candidate_regions],
            "recommended_positions": unique_positions,
            "status": "AVAILABLE",
            "reviewed": False,  # Mandatory human review invariant
            "provenance_basis": "CURATED_UNIPROT_RECORD",
            "biological_validation": "NOT_AVAILABLE",
            "epistemic_disclaimer": (
                "UniProt feature records describe native biological properties. Translating these "
                "features into translation pause recommendations is an exploratory computational "
                "heuristic; biological validation: NOT_AVAILABLE."
            ),
        }

        return sidecar_data


def asdict_cand(cr: CandidateRegion) -> Dict[str, Any]:
    return {
        "start_1based": cr.start_1based,
        "end_1based": cr.end_1based,
        "length": cr.length,
        "region_type": cr.region_type,
        "reason": cr.reason,
        "heuristic_basis": cr.heuristic_basis,
    }


def fetch_uniprot_sidecar(
    accession: str,
    target_protein_sequence: Optional[str] = None,
    cache_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Convenience top-level functional interface for UniProt sidecar generation."""
    fetcher = UniProtSidecarFetcher(cache_dir=cache_dir)
    return fetcher.fetch_and_format_sidecar(
        accession=accession,
        target_protein_sequence=target_protein_sequence,
    )
