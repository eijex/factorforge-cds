"""Shared sequence metrics for CDS evaluation and validation."""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from factorforge.engines.profile.utils import get_data_path

# Homopolymer thresholds — two distinct concerns, intentionally different values.
# Expression stability: AT-rich runs ≥6 nt can resemble instability elements (ARE).
# Synthesis/manufacturing: runs ≥8 nt are flagged by gene synthesis vendors as difficult.
HOMOPOLYMER_EXPRESSION_WARN_NT = 6
HOMOPOLYMER_SYNTHESIS_WARN_NT = 8


STANDARD_GENETIC_CODE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

STOP_CODONS = {codon for codon, aa in STANDARD_GENETIC_CODE.items() if aa == "*"}
VALID_BASES = set("ATGC")


@dataclass(frozen=True)
class CodonUsageTable:
    codon_to_aa: dict[str, str]
    codon_weights: dict[str, float]
    best_codon_for_aa: dict[str, str]
    source: str | None = None


def _default_codon_table_path() -> Path:
    # Job 168 / v3.3.0 (_analysis/025): production default switched from the
    # legacy Kazusa/SGN-derived table to the NbeV1.1 LAB-strain high-confidence
    # derived table. See data/reference/active_codon_reference.json.
    return get_data_path() / "profiles" / "nbev11_cds_hc_derived_codons.json"


def load_codon_usage_table(path: Path | None = None) -> CodonUsageTable:
    table_path = path or _default_codon_table_path()
    raw = json.loads(table_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Codon table JSON must be an object")

    codons = raw.get("codons")
    if not isinstance(codons, dict):
        raise ValueError("Codon table JSON missing 'codons'")

    codon_to_aa: dict[str, str] = {}
    codon_freq: dict[str, float] = {}
    for codon, entry in codons.items():
        if not isinstance(entry, dict) or not isinstance(codon, str):
            continue
        aa = entry.get("aa")
        freq = entry.get("frequency")
        if not isinstance(aa, str) or not isinstance(freq, (int, float)):
            continue
        codon_to_aa[codon] = aa
        codon_freq[codon] = float(freq)

    codon_weights = _build_codon_weights(codon_to_aa, codon_freq)
    best_codon_for_aa = _best_codon_map(codon_to_aa, codon_freq)
    source = raw.get("source") if isinstance(raw.get("source"), str) else None

    return CodonUsageTable(
        codon_to_aa=codon_to_aa,
        codon_weights=codon_weights,
        best_codon_for_aa=best_codon_for_aa,
        source=source,
    )


def _build_codon_weights(
    codon_to_aa: dict[str, str],
    codon_freq: dict[str, float],
) -> dict[str, float]:
    by_aa: dict[str, list[float]] = {}
    for codon, aa in codon_to_aa.items():
        if aa == "*":
            continue
        by_aa.setdefault(aa, []).append(codon_freq.get(codon, 0.0))

    weights: dict[str, float] = {}
    for codon, aa in codon_to_aa.items():
        if aa == "*":
            continue
        max_freq = max(by_aa.get(aa, [0.0]))
        freq = codon_freq.get(codon, 0.0)
        weights[codon] = freq / max_freq if max_freq > 0 else 0.0
    return weights


def _best_codon_map(
    codon_to_aa: dict[str, str],
    codon_freq: dict[str, float],
) -> dict[str, str]:
    best: dict[str, tuple[str, float]] = {}
    for codon, aa in codon_to_aa.items():
        if aa == "*":
            continue
        current = best.get(aa)
        freq = codon_freq.get(codon, 0.0)
        if current is None or freq > current[1]:
            best[aa] = (codon, freq)
    return {aa: codon for aa, (codon, _) in best.items()}


def _normalize_dna(sequence: str) -> str:
    return "".join(sequence.upper().replace("U", "T").split())


def _codons(sequence: str, include_partial: bool = False) -> list[str]:
    seq = _normalize_dna(sequence)
    end = len(seq) if include_partial else len(seq) - len(seq) % 3
    return [seq[index : index + 3] for index in range(0, end, 3)]


def calculate_gc(sequence: str) -> float:
    """Calculate GC content as a percentage in the range 0-100."""
    seq = _normalize_dna(sequence)
    if not seq:
        return 0.0
    return ((seq.count("G") + seq.count("C")) / len(seq)) * 100.0


def calculate_gc_windows(
    sequence: str,
    window_size: int = 60,
    step: int = 30,
) -> list[dict[str, float | int]]:
    """Calculate sliding-window GC percentages."""
    if window_size <= 0:
        raise ValueError("window_size must be > 0")
    if step <= 0:
        raise ValueError("step must be > 0")

    seq = _normalize_dna(sequence)
    if not seq:
        return []
    if len(seq) <= window_size:
        return [{"start": 0, "end": len(seq), "gc": calculate_gc(seq)}]

    windows: list[dict[str, float | int]] = []
    for start in range(0, len(seq) - window_size + 1, step):
        end = start + window_size
        windows.append({"start": start, "end": end, "gc": calculate_gc(seq[start:end])})
    if windows and int(windows[-1]["end"]) < len(seq):
        start = len(seq) - window_size
        if start != int(windows[-1]["start"]):
            windows.append({"start": start, "end": len(seq), "gc": calculate_gc(seq[start:])})
    return windows


def calculate_first_region_gc(
    sequence: str,
    region_sizes: list[int] | None = None,
) -> dict[str, float]:
    """Calculate GC for configured 5-prime regions."""
    seq = _normalize_dna(sequence)
    sizes = region_sizes or [30, 60, 90]
    result: dict[str, float] = {}
    for size in sizes:
        if size <= 0:
            raise ValueError("region sizes must be > 0")
        region = seq[: min(size, len(seq))]
        result[f"first_{size}nt_gc"] = calculate_gc(region)
    return result


def translate_dna(sequence: str) -> str:
    """Translate DNA to amino acids, using X for invalid codons."""
    translated: list[str] = []
    for codon in _codons(sequence):
        translated.append(STANDARD_GENETIC_CODE.get(codon, "X"))
    return "".join(translated)


def amino_acid_identity(input_protein: str, dna_sequence: str) -> float:
    """Return exact-position amino acid identity from translated DNA."""
    expected = "".join(input_protein.upper().split())
    observed = translate_dna(dna_sequence)
    if observed.endswith("*") and not expected.endswith("*"):
        observed = observed[:-1]
    if not expected:
        return 0.0
    matches = sum(1 for exp, obs in zip(expected, observed) if exp == obs)
    return matches / len(expected) if len(observed) == len(expected) else matches / len(expected)


def count_internal_stops(dna_sequence: str) -> int:
    """Count stop codons before the final codon."""
    codons = _codons(dna_sequence)
    return sum(1 for codon in codons[:-1] if codon in STOP_CODONS)


def calculate_cai(sequence: str, codon_weights: dict[str, float]) -> float:
    """Calculate CAI as a geometric mean of supplied codon weights."""
    seq = _normalize_dna(sequence)
    if not seq or len(seq) % 3 != 0:
        return 0.0

    log_sum = 0.0
    count = 0
    for codon in _codons(seq):
        if codon in STOP_CODONS:
            continue
        weight = codon_weights.get(codon)
        if weight is None or weight <= 0:
            return 0.0
        log_sum += math.log(weight)
        count += 1
    return math.exp(log_sum / count) if count else 0.0


def calculate_dinucleotide_score(
    sequence: str,
    cpg_weight: float = 0.0,
    tpa_weight: float = 1.0,
) -> float:
    """Score dinucleotide avoidance.

    Plant default: CpG inactive (cpg_weight=0.0), only TpA is penalized.
    Mammalian opt-in: set cpg_weight=1.0 and tpa_weight=1.0 to penalize both.
    """
    from factorforge.engines.profile.utils import calculate_dinucleotide_ratio

    if len(sequence) < 6:
        return 1.0

    total_weight = cpg_weight + tpa_weight
    if total_weight == 0:
        return 1.0

    cpg_ratio = calculate_dinucleotide_ratio(sequence, "CG")
    tpa_ratio = calculate_dinucleotide_ratio(sequence, "TA")
    cpg_score = max(0.0, 1.0 - cpg_ratio / 2.0)
    tpa_score = max(0.0, 1.0 - tpa_ratio / 2.0)
    return (cpg_weight * cpg_score + tpa_weight * tpa_score) / total_weight


def codon_usage_profile(sequence: str) -> dict[str, dict[str, float | int | str]]:
    """Return codon counts and frequencies for a DNA sequence."""
    codons = _codons(sequence)
    counts = Counter(codons)
    total = sum(counts.values())
    profile: dict[str, dict[str, float | int | str]] = {}
    for codon in sorted(counts):
        profile[codon] = {
            "count": counts[codon],
            "frequency": counts[codon] / total if total else 0.0,
            "aa": STANDARD_GENETIC_CODE.get(codon, "X"),
        }
    return profile


def detect_homopolymers(
    sequence: str,
    max_run: int = HOMOPOLYMER_EXPRESSION_WARN_NT,
) -> list[dict[str, Any]]:
    """Detect homopolymer runs for expression stability evaluation.

    Uses HOMOPOLYMER_EXPRESSION_WARN_NT (default 6 nt) — AT-rich runs of this
    length can resemble AU-rich instability elements (ARE) and affect mRNA
    stability in plant expression systems.

    For synthesis/manufacturing risk, see RuleEngine.scan_homopolymers()
    which uses HOMOPOLYMER_SYNTHESIS_WARN_NT (8 nt).
    """
    if max_run <= 1:
        raise ValueError("max_run must be > 1")

    seq = _normalize_dna(sequence)
    findings: list[dict[str, Any]] = []
    if not seq:
        return findings

    run_base = seq[0]
    run_start = 0
    for index, base in enumerate(seq[1:], start=1):
        if base == run_base:
            continue
        run_length = index - run_start
        if run_length >= max_run:
            findings.append({
                "start": run_start,
                "end": index,
                "base": run_base,
                "length": run_length,
                "context": "expression_stability",
                "threshold_nt": max_run,
            })
        run_base = base
        run_start = index

    run_length = len(seq) - run_start
    if run_length >= max_run:
        findings.append({
            "start": run_start,
            "end": len(seq),
            "base": run_base,
            "length": run_length,
            "context": "expression_stability",
            "threshold_nt": max_run,
        })
    return findings


def detect_repeats(sequence: str) -> list[dict[str, Any]]:
    """Detect simple tandem repeats with motif length 2-6 repeated at least 3 times."""
    seq = _normalize_dna(sequence)
    findings: list[dict[str, Any]] = []
    occupied: set[tuple[int, int]] = set()

    for motif_len in range(2, 7):
        index = 0
        while index <= len(seq) - motif_len * 3:
            motif = seq[index : index + motif_len]
            repeats = 1
            cursor = index + motif_len
            while seq[cursor : cursor + motif_len] == motif:
                repeats += 1
                cursor += motif_len
            if repeats >= 3:
                span = (index, cursor)
                if span not in occupied:
                    findings.append(
                        {
                            "start": index,
                            "end": cursor,
                            "motif": motif,
                            "repeat_count": repeats,
                        }
                    )
                    occupied.add(span)
                index = cursor
            else:
                index += 1
    return findings


def detect_forbidden_motifs(sequence: str, motifs: list[str]) -> list[dict[str, Any]]:
    """Find all exact forbidden motif occurrences."""
    seq = _normalize_dna(sequence)
    findings: list[dict[str, Any]] = []
    for motif in motifs:
        normalized = _normalize_dna(motif)
        if not normalized:
            continue
        start = seq.find(normalized)
        while start != -1:
            findings.append({"start": start, "end": start + len(normalized), "motif": normalized})
            start = seq.find(normalized, start + 1)
    return findings


def detect_invalid_codons(sequence: str) -> list[dict[str, Any]]:
    """Detect invalid, partial, or non-ATGC codons."""
    seq = _normalize_dna(sequence)
    findings: list[dict[str, Any]] = []
    for index, codon in enumerate(_codons(seq, include_partial=True)):
        start = index * 3
        if len(codon) != 3:
            findings.append({"start": start, "end": len(seq), "codon": codon, "reason": "partial"})
        elif set(codon) - VALID_BASES:
            findings.append({"start": start, "end": start + 3, "codon": codon, "reason": "invalid_base"})
        elif codon not in STANDARD_GENETIC_CODE:
            findings.append({"start": start, "end": start + 3, "codon": codon, "reason": "unknown"})
    return findings
