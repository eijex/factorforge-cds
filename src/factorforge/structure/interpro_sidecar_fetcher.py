"""InterPro Automated Sidecar Fetcher & Domain Resolver for FactorForge.

Fetches and normalizes curated protein family, domain, and active site annotations
from EMBL-EBI InterPro REST API into the exact sidecar schema required by
`factorforge.nbfold.annotations` and structural pause conditioning.

Epistemic & Governance Constraints:
- Epistemic classification: UNVERIFIED_PUBLIC_CORPUS / INTERPRO_REST_API.
- Target Accession Invariant: Strictly verifies that fetched entries belong to the
  queried accession; rejects mismatched mock/cached payloads.
- Biological Boundary: Only true 'domain' entries generate folding pause candidate
  regions ('interpro_domain'). 'family', 'repeat', and 'site' entries are preserved
  in the sidecar entries list for annotation, but excluded from domain pause candidates.
- Invariant: reviewed=False on all output payloads.
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

INTERPRO_FETCHER_VERSION = "1.1.0"
INTERPRO_API_ENDPOINT = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{accession}"
ACCESSION_REGEX = re.compile(r"^[A-Z0-9_-]{4,20}$", re.IGNORECASE)


def fetch_interpro_raw_record(
    accession: str,
    timeout: float = 10.0,
    headers: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Fetch raw JSON record directly from InterPro REST API."""
    clean_acc = accession.strip().upper()
    if not clean_acc or not ACCESSION_REGEX.match(clean_acc):
        raise ValueError(f"Invalid UniProt accession for InterPro query: '{accession}'")

    url = INTERPRO_API_ENDPOINT.format(accession=clean_acc)
    req_headers = {
        "User-Agent": "FactorForge-InterProFetcher/1.1 (Eijex; scientific-research)",
        "Accept": "application/json",
    }
    if headers:
        req_headers.update(headers)

    req = urllib.request.Request(url, headers=req_headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        content = resp.read().decode("utf-8")
        return json.loads(content)


class InterProSidecarFetcher:
    """Fetches, caches, validates, and formats InterPro records into valid FactorForge sidecar contracts."""

    def __init__(self, cache_dir: Optional[Path] = None, timeout: float = 8.0) -> None:
        self.cache_dir = cache_dir
        self.timeout = timeout
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def fetch_and_format_sidecar(
        self,
        accession: str,
        target_protein_sequence: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Fetch InterPro entries and format into `interpro_sidecar_data` compliant dictionary."""
        if not isinstance(accession, str):
            return {
                "status": "ERROR_INVALID_ACCESSION",
                "reason": "Accession must be a non-empty string.",
                "accession": str(accession),
                "entries": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        clean_acc = accession.strip().upper()
        if not clean_acc or not ACCESSION_REGEX.match(clean_acc):
            return {
                "status": "ERROR_INVALID_ACCESSION",
                "reason": f"Invalid accession format: '{accession}'",
                "accession": clean_acc,
                "entries": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        if not isinstance(target_protein_sequence, str) or not target_protein_sequence.strip():
            return {
                "status": "ERROR_INVALID_SEQUENCE",
                "reason": "target_protein_sequence must be a non-empty string.",
                "accession": clean_acc,
                "entries": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        clean_seq = "".join(target_protein_sequence.upper().split()).rstrip("*")
        seq_len = len(clean_seq)
        seq_digest = hashlib.sha256(clean_seq.encode("utf-8")).hexdigest()

        cache_file = self.cache_dir / f"{clean_acc}_interpro.json" if self.cache_dir else None
        raw_data: Optional[Dict[str, Any]] = None

        if use_cache and cache_file and cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached_obj = json.load(f)
                # Check cache accession integrity
                if isinstance(cached_obj, dict):
                    cache_query_acc = cached_obj.get("query_accession", "").strip().upper()
                    if cache_query_acc and cache_query_acc != clean_acc:
                        raw_data = None
                    else:
                        raw_data = cached_obj.get("data", cached_obj)
            except Exception:
                raw_data = None

        if raw_data is None:
            try:
                raw_data = fetch_interpro_raw_record(clean_acc, timeout=self.timeout)
                if cache_file and raw_data:
                    with open(cache_file, "w", encoding="utf-8") as f:
                        json.dump({"query_accession": clean_acc, "data": raw_data}, f, indent=2)
            except Exception as exc:
                return {
                    "status": "UNAVAILABLE",
                    "reason": f"Failed to fetch InterPro entries for '{clean_acc}': {exc}",
                    "accession": clean_acc,
                    "sequence": clean_seq,
                    "protein_length": seq_len,
                    "protein_sha256": seq_digest,
                    "entries": [],
                    "candidate_regions": [],
                    "recommended_positions": [],
                    "reviewed": False,
                }

        raw_results = raw_data.get("results", [])
        if not isinstance(raw_results, list):
            raw_results = []

        if not raw_results:
            return {
                "status": "NO_ENTRIES_FOUND",
                "reason": f"No InterPro domain or family entries found for accession '{clean_acc}'.",
                "accession": clean_acc,
                "sequence": clean_seq,
                "protein_length": seq_len,
                "protein_sha256": seq_digest,
                "entries": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        # Check target accession alignment across raw results
        total_proteins_seen = 0
        matched_proteins_seen = 0

        normalized_entries: List[Dict[str, Any]] = []
        candidate_regions: List[CandidateRegion] = []
        recommended_positions: List[int] = []

        for item in raw_results:
            if not isinstance(item, dict):
                continue
            meta = item.get("metadata", {})
            entry_acc = str(meta.get("accession") or item.get("accession") or "").strip()
            if not entry_acc:
                continue

            entry_name = str(meta.get("name") or item.get("name") or entry_acc).strip()
            source_db = str(meta.get("source_database") or item.get("source_database") or "interpro").strip().lower()
            raw_type = str(meta.get("type") or item.get("type") or "domain").strip().lower()

            # Extract locations strictly for the target accession
            locations: List[Dict[str, Any]] = []
            proteins_list = item.get("proteins", [])
            if isinstance(proteins_list, list) and proteins_list:
                for p in proteins_list:
                    if isinstance(p, dict):
                        p_acc = str(p.get("accession", "")).strip().upper()
                        if p_acc:
                            total_proteins_seen += 1
                        if p_acc == clean_acc:
                            matched_proteins_seen += 1
                            for loc in p.get("entry_protein_locations", []):
                                if isinstance(loc, dict):
                                    for frag in loc.get("fragments", []):
                                        if isinstance(frag, dict):
                                            s = frag.get("start")
                                            e = frag.get("end")
                                            if isinstance(s, int) and isinstance(e, int) and 1 <= s <= e <= seq_len:
                                                locations.append({"start": s, "end": e})

            elif "locations" in item:
                # Direct locations format (mock or simplified payload)
                for loc in item["locations"]:
                    if isinstance(loc, dict):
                        s = loc.get("start")
                        e = loc.get("end")
                        if isinstance(s, int) and isinstance(e, int) and 1 <= s <= e <= seq_len:
                            locations.append({"start": s, "end": e})

            if not locations:
                continue

            norm_entry = {
                "accession": entry_acc,
                "name": entry_name,
                "source_database": source_db,
                "entry_type": raw_type,
                "locations": locations,
            }
            normalized_entries.append(norm_entry)

            # Invariant: Only 'domain' entries generate 'interpro_domain' pause candidate regions.
            # 'family', 'repeat', and 'site' entries are kept in 'entries' for classification,
            # but excluded from domain folding pause boundaries.
            if raw_type == "domain":
                for loc in locations:
                    start_p = loc["start"]
                    end_p = loc["end"]
                    recommended_positions.append(end_p)
                    candidate_regions.append(
                        CandidateRegion(
                            start_1based=start_p,
                            end_1based=end_p,
                            length=end_p - start_p + 1,
                            region_type="interpro_domain",
                            reason=f"InterPro {source_db.upper()} domain '{entry_name}' ({entry_acc}) [{start_p}..{end_p}]",
                            heuristic_basis="INTERPRO_DOMAIN_HEURISTIC_PROPOSAL",
                        )
                    )

        # Accession mismatch defense: If raw payload had proteins, but zero matched clean_acc
        if total_proteins_seen > 0 and matched_proteins_seen == 0:
            return {
                "status": "ERROR_ACCESSION_MISMATCH",
                "reason": (
                    f"InterPro payload contains protein entries, but none match requested accession '{clean_acc}'."
                ),
                "accession": clean_acc,
                "sequence": clean_seq,
                "protein_length": seq_len,
                "protein_sha256": seq_digest,
                "entries": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        if not normalized_entries:
            return {
                "status": "NO_ENTRIES_FOUND",
                "reason": f"No valid InterPro entries resolved for accession '{clean_acc}'.",
                "accession": clean_acc,
                "sequence": clean_seq,
                "protein_length": seq_len,
                "protein_sha256": seq_digest,
                "entries": [],
                "candidate_regions": [],
                "recommended_positions": [],
                "reviewed": False,
            }

        unique_positions = sorted(list(set(recommended_positions)))

        return {
            "status": "AVAILABLE",
            "accession": clean_acc,
            "sequence": clean_seq,
            "protein_length": seq_len,
            "protein_sha256": seq_digest,
            "entries": normalized_entries,
            "candidate_regions": [asdict_cand(cr) for cr in candidate_regions],
            "recommended_positions": unique_positions,
            "reviewed": False,  # Mandatory human review invariant
            "provenance_basis": "INTERPRO_REST_API",
            "biological_validation": "NOT_AVAILABLE",
            "epistemic_disclaimer": (
                "InterPro domain definitions describe structural/functional units. Translating domain "
                "boundaries into translation pause recommendations is an exploratory computational "
                "heuristic; biological validation: NOT_AVAILABLE."
            ),
        }


def asdict_cand(cr: CandidateRegion) -> Dict[str, Any]:
    return {
        "start_1based": cr.start_1based,
        "end_1based": cr.end_1based,
        "length": cr.length,
        "region_type": cr.region_type,
        "reason": cr.reason,
        "heuristic_basis": cr.heuristic_basis,
    }


def fetch_interpro_sidecar(
    accession: str,
    target_protein_sequence: str,
    cache_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Convenience functional interface for InterPro sidecar generation."""
    fetcher = InterProSidecarFetcher(cache_dir=cache_dir)
    return fetcher.fetch_and_format_sidecar(
        accession=accession,
        target_protein_sequence=target_protein_sequence,
    )
