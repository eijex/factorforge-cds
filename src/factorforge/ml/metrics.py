"""Shared sequence metrics for v3-alpha evaluation and validation."""

from __future__ import annotations

import math
from collections import Counter
from typing import Any


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


def detect_homopolymers(sequence: str, max_run: int = 6) -> list[dict[str, Any]]:
    """Detect runs whose length is greater than or equal to max_run."""
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
            findings.append(
                {"start": run_start, "end": index, "base": run_base, "length": run_length}
            )
        run_base = base
        run_start = index

    run_length = len(seq) - run_start
    if run_length >= max_run:
        findings.append(
            {"start": run_start, "end": len(seq), "base": run_base, "length": run_length}
        )
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

