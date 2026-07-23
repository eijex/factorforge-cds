"""Deterministic CDS design-review parsing, criteria, and decision helpers."""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import re
from typing import Any, Literal

from factorforge.analysis.metrics import (
    STOP_CODONS,
    calculate_cai,
    calculate_gc,
    calculate_gc_windows,
    detect_forbidden_motifs,
    translate_dna,
)

ConstraintMode = Literal["required", "preferred", "ignored"]
AutomatedDecision = Literal["PASS", "CONDITIONAL_PASS", "FAIL"]

MULTI_FASTA_ERROR = "Multiple FASTA records detected. Upload one sequence at a time."
VALID_PROTEIN = frozenset("ACDEFGHIKLMNPQRSTVWY*")
VALID_DNA = frozenset("ACGT")
TYPE_IIS_SITES = {
    "BsaI": ("GGTCTC", "GAGACC"),
    "BsmBI/Esp3I": ("CGTCTC", "GAGACG"),
    "SapI": ("GCTCTTC", "GAAGAGC"),
}
VALID_MODES = frozenset({"required", "preferred", "ignored"})
VALID_DISPOSITIONS = frozenset(
    {"accept", "accept_with_exception", "return_for_redesign", "reject"}
)


def _compact_sequence(value: str) -> str:
    return re.sub(r"\s+", "", value).upper()


def parse_sequence_input(raw: str) -> dict[str, Any]:
    """Parse one plain sequence or one FASTA record and classify it explicitly."""
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError("Sequence is required.")

    lines = raw.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    nonempty = [line.strip() for line in lines if line.strip()]
    headers = [index for index, line in enumerate(nonempty) if line.startswith(">")]
    if len(headers) > 1:
        raise ValueError(MULTI_FASTA_ERROR)
    if headers and headers[0] != 0:
        raise ValueError("FASTA header must be the first non-empty line.")

    header = nonempty[0][1:].strip() if headers else None
    sequence_lines = nonempty[1:] if headers else nonempty
    sequence = _compact_sequence("".join(sequence_lines))
    if not sequence:
        raise ValueError("Sequence is required.")

    characters = set(sequence)
    if characters <= VALID_DNA:
        return _build_cds_context(sequence, header)

    invalid = characters - VALID_PROTEIN
    if invalid:
        raise ValueError(f"Invalid sequence characters: {', '.join(sorted(invalid))}.")
    if "*" in sequence[:-1]:
        raise ValueError("Protein input contains an internal stop marker.")

    protein = sequence.rstrip("*")
    if not protein:
        raise ValueError("Protein sequence must contain at least one amino acid.")
    return {
        "input_type": "protein",
        "normalized_sequence": sequence,
        "optimization_sequence": protein,
        "fasta_header": header,
        "generation_allowed": True,
        "message": (
            "Protein sequence detected. FactorForge will reverse translate it into "
            "a host-adapted coding sequence."
        ),
        "summary": {
            "protein_length_aa": len(protein),
            "terminal_stop_marker": sequence.endswith("*"),
        },
        "errors": [],
        "warnings": [],
    }


def _build_cds_context(sequence: str, header: str | None) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    frame_valid = len(sequence) % 3 == 0
    if not frame_valid:
        errors.append("CDS length must be divisible by 3.")

    codons = [sequence[index : index + 3] for index in range(0, len(sequence), 3)]
    complete_codons = [codon for codon in codons if len(codon) == 3]
    start_present = bool(complete_codons) and complete_codons[0] == "ATG"
    terminal_stop_present = bool(complete_codons) and complete_codons[-1] in STOP_CODONS
    internal_stops = [
        index for index, codon in enumerate(complete_codons[:-1]) if codon in STOP_CODONS
    ]
    if internal_stops:
        errors.append("CDS contains internal stop codon(s).")
    if not start_present:
        warnings.append("Start codon is absent.")
    if not terminal_stop_present:
        warnings.append("Terminal stop codon is absent.")

    translated = translate_dna(sequence) if frame_valid else ""
    protein = translated[:-1] if terminal_stop_present else translated
    return {
        "input_type": "cds",
        "normalized_sequence": sequence,
        "optimization_sequence": protein,
        "fasta_header": header,
        "generation_allowed": not errors,
        "message": (
            "DNA coding sequence detected. FactorForge will preserve the translated "
            "protein while redesigning synonymous codons for the selected host and constraints."
        ),
        "summary": {
            "cds_length_bp": len(sequence),
            "codon_count": len(complete_codons),
            "protein_length_aa": len(protein),
            "gc_percent": round(calculate_gc(sequence), 1),
            "frame_valid": frame_valid,
            "start_codon_present": start_present,
            "terminal_stop_present": terminal_stop_present,
            "terminal_stop_codon": complete_codons[-1] if terminal_stop_present else None,
            "internal_stop_count": len(internal_stops),
        },
        "errors": errors,
        "warnings": warnings,
    }


def restore_cds_stop_policy(candidate: str, input_context: dict[str, Any]) -> str:
    """Restore the original CDS terminal-stop policy to a generated candidate."""
    sequence = _compact_sequence(candidate)
    if input_context["input_type"] != "cds":
        return sequence
    stop = input_context["summary"]["terminal_stop_codon"]
    translated = translate_dna(sequence)
    if translated.endswith("*"):
        sequence = sequence[:-3]
    return sequence + (stop or "")


def assert_pathway_invariants(
    original_cds: str,
    optimized_cds: str,
    input_context: dict[str, Any],
) -> None:
    """Raise visibly when CDS length or translated-protein invariants regress."""
    if input_context["input_type"] != "cds":
        return
    if len(original_cds) != len(optimized_cds):
        raise ValueError(
            "CDS invariant failed: optimized nucleotide length differs from the original."
        )
    original_protein = translate_dna(original_cds).rstrip("*")
    optimized_protein = translate_dna(optimized_cds).rstrip("*")
    if original_protein != optimized_protein:
        raise ValueError("CDS invariant failed: translated protein was not preserved.")


def default_acceptance_criteria(
    gc_min: float = 40.0,
    gc_max: float = 55.0,
) -> dict[str, dict[str, Any]]:
    return {
        "cai": {"mode": "preferred", "minimum": 0.8},
        "overall_gc": {"mode": "preferred", "minimum": gc_min, "maximum": gc_max},
        "local_gc": {
            "mode": "preferred",
            "minimum": 30.0,
            "maximum": 70.0,
            "window_size": 60,
        },
        "type_iis": {
            "mode": "required",
            "enzymes": ["BsaI", "BsmBI/Esp3I", "SapI"],
            "custom_sites": [],
        },
        "repeats": {"mode": "preferred", "minimum_length": 18, "maximum_count": 0},
        "homopolymers": {"mode": "preferred", "maximum_length": 8},
        "forbidden_motifs": {
            "mode": "preferred",
            "motifs": ["AATAAA", "GTAAGT", "ATTTA"],
        },
    }


def normalize_acceptance_criteria(
    value: dict[str, Any] | None,
    *,
    gc_min: float,
    gc_max: float,
) -> dict[str, dict[str, Any]]:
    criteria = default_acceptance_criteria(gc_min, gc_max)
    if value is None:
        return criteria
    if not isinstance(value, dict):
        raise ValueError("acceptance_criteria must be an object.")
    unknown = set(value) - set(criteria)
    if unknown:
        raise ValueError(f"Unknown acceptance criteria: {', '.join(sorted(unknown))}.")
    for name, override in value.items():
        if not isinstance(override, dict):
            raise ValueError(f"acceptance_criteria.{name} must be an object.")
        criteria[name].update(deepcopy(override))
        mode = criteria[name].get("mode")
        if mode not in VALID_MODES:
            raise ValueError(
                f"acceptance_criteria.{name}.mode must be required, preferred, or ignored."
            )
    _validate_criteria_ranges(criteria)
    return criteria


def _validate_criteria_ranges(criteria: dict[str, dict[str, Any]]) -> None:
    for name in ("overall_gc", "local_gc"):
        minimum = float(criteria[name]["minimum"])
        maximum = float(criteria[name]["maximum"])
        if not 0 <= minimum <= maximum <= 100:
            raise ValueError(f"acceptance_criteria.{name} range must be within 0-100.")
    cai_minimum = float(criteria["cai"]["minimum"])
    if not 0 <= cai_minimum <= 1:
        raise ValueError("acceptance_criteria.cai.minimum must be within 0-1.")
    if int(criteria["local_gc"]["window_size"]) <= 0:
        raise ValueError("acceptance_criteria.local_gc.window_size must be positive.")
    if int(criteria["homopolymers"]["maximum_length"]) < 2:
        raise ValueError("acceptance_criteria.homopolymers.maximum_length must be at least 2.")


def _count_direct_repeats(sequence: str, minimum_length: int) -> int:
    if len(sequence) < minimum_length * 2:
        return 0
    seen: set[str] = set()
    repeated: set[str] = set()
    for index in range(0, len(sequence) - minimum_length + 1):
        motif = sequence[index : index + minimum_length]
        if motif in seen:
            repeated.add(motif)
        seen.add(motif)
    return len(repeated)


def _longest_homopolymer(sequence: str) -> int:
    return max((len(match.group(0)) for match in re.finditer(r"([ACGT])\1*", sequence)), default=0)


def _type_iis_findings(
    sequence: str,
    criterion: dict[str, Any],
) -> list[dict[str, Any]]:
    findings: list[dict[str, Any]] = []
    enabled = set(criterion.get("enzymes", []))
    for enzyme, patterns in TYPE_IIS_SITES.items():
        if enzyme not in enabled:
            continue
        for pattern in patterns:
            start = sequence.find(pattern)
            while start >= 0:
                findings.append({"enzyme": enzyme, "site": pattern, "start": start})
                start = sequence.find(pattern, start + 1)
    for site in criterion.get("custom_sites", []):
        if isinstance(site, dict):
            name = str(site.get("name") or "Custom")
            pattern = _compact_sequence(str(site.get("sequence") or ""))
        else:
            name = "Custom"
            pattern = _compact_sequence(str(site))
        if not pattern or set(pattern) - VALID_DNA:
            continue
        start = sequence.find(pattern)
        while start >= 0:
            findings.append({"enzyme": name, "site": pattern, "start": start})
            start = sequence.find(pattern, start + 1)
    return findings


def evaluate_candidate(
    sequence: str,
    *,
    cai: float,
    criteria: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Evaluate one CDS against a normalized criteria snapshot."""
    seq = _compact_sequence(sequence)
    overall_gc = calculate_gc(seq)
    local_rule = criteria["local_gc"]
    windows = calculate_gc_windows(
        seq,
        window_size=int(local_rule["window_size"]),
        step=max(1, int(local_rule["window_size"]) // 2),
    )
    local_values = [float(window["gc"]) for window in windows] or [overall_gc]
    type_iis = _type_iis_findings(seq, criteria["type_iis"])
    repeat_count = _count_direct_repeats(seq, int(criteria["repeats"]["minimum_length"]))
    longest_homopolymer = _longest_homopolymer(seq)
    motifs = detect_forbidden_motifs(seq, list(criteria["forbidden_motifs"].get("motifs", [])))

    observations = {
        "cai": (float(cai), float(criteria["cai"]["minimum"]), float(cai) >= float(criteria["cai"]["minimum"])),
        "overall_gc": (
            round(overall_gc, 1),
            f"{float(criteria['overall_gc']['minimum']):g}-{float(criteria['overall_gc']['maximum']):g}%",
            float(criteria["overall_gc"]["minimum"]) <= overall_gc <= float(criteria["overall_gc"]["maximum"]),
        ),
        "local_gc": (
            f"{min(local_values):.1f}-{max(local_values):.1f}%",
            f"{float(local_rule['minimum']):g}-{float(local_rule['maximum']):g}%",
            min(local_values) >= float(local_rule["minimum"]) and max(local_values) <= float(local_rule["maximum"]),
        ),
        "type_iis": (len(type_iis), 0, not type_iis),
        "repeats": (
            repeat_count,
            int(criteria["repeats"]["maximum_count"]),
            repeat_count <= int(criteria["repeats"]["maximum_count"]),
        ),
        "homopolymers": (
            longest_homopolymer,
            int(criteria["homopolymers"]["maximum_length"]),
            longest_homopolymer <= int(criteria["homopolymers"]["maximum_length"]),
        ),
        "forbidden_motifs": (len(motifs), 0, not motifs),
    }

    rows: list[dict[str, Any]] = []
    required_failures = 0
    preferred_warnings = 0
    for name, (observed, threshold, passed) in observations.items():
        mode: ConstraintMode = criteria[name]["mode"]
        if mode == "ignored":
            result = "IGNORED"
        elif passed:
            result = "PASS"
        elif mode == "required":
            result = "FAIL"
            required_failures += 1
        else:
            result = "WARN"
            preferred_warnings += 1
        rows.append(
            {
                "criterion": name,
                "mode": mode,
                "observed": observed,
                "threshold": threshold,
                "result": result,
            }
        )

    if required_failures:
        decision: AutomatedDecision = "FAIL"
    elif preferred_warnings:
        decision = "CONDITIONAL_PASS"
    else:
        decision = "PASS"
    failed_names = [
        row["criterion"] for row in rows if row["result"] in {"FAIL", "WARN"}
    ]
    explanation = (
        f"{', '.join(failed_names)} require review."
        if failed_names
        else "All active acceptance criteria passed."
    )
    return {
        "automated_decision": decision,
        "required_failure_count": required_failures,
        "preferred_warning_count": preferred_warnings,
        "explanation": explanation,
        "criteria": rows,
        "details": {
            "type_iis_sites": type_iis,
            "forbidden_motifs": motifs,
            "local_gc_min": round(min(local_values), 1),
            "local_gc_max": round(max(local_values), 1),
            "repeat_count": repeat_count,
            "longest_homopolymer": longest_homopolymer,
        },
    }


def apply_reviewer_disposition(
    automated_decision: AutomatedDecision,
    value: dict[str, Any] | None,
) -> dict[str, Any] | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("reviewer_disposition must be an object.")
    disposition = str(value.get("disposition") or "").strip().lower()
    if disposition not in VALID_DISPOSITIONS:
        raise ValueError(
            "reviewer_disposition.disposition must be accept, accept_with_exception, "
            "return_for_redesign, or reject."
        )
    reason = str(value.get("reason") or "").strip()
    accepting_failure = automated_decision == "FAIL" and disposition in {
        "accept",
        "accept_with_exception",
    }
    if accepting_failure and not reason:
        raise ValueError("A written reason is required to accept an automated FAIL.")
    final_state = {
        "accept": "ACCEPTED",
        "accept_with_exception": "ACCEPTED_WITH_EXCEPTION",
        "return_for_redesign": "RETURNED_FOR_REDESIGN",
        "reject": "REJECTED",
    }[disposition]
    if accepting_failure:
        final_state = "MANUALLY_ACCEPTED"
    return {
        "disposition": disposition,
        "reason": reason or None,
        "timestamp": value.get("timestamp")
        or datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "final_state": final_state,
        "automated_decision": automated_decision,
    }


def build_result_identifier(input_sequence: str, criteria: dict[str, Any]) -> str:
    payload = f"{input_sequence}|{repr(criteria)}".encode("utf-8")
    return "ff-" + hashlib.sha256(payload).hexdigest()[:16]


def original_cai(sequence: str, codon_weights: dict[str, float]) -> float | None:
    if not sequence or len(sequence) % 3:
        return None
    return calculate_cai(sequence, codon_weights)
