"""Equivalence-preserving FactorForge Default SOP refinement.

Applies:
1. Optional 5' signal-peptide prefix preservation (first 60 nt / 20 codons)
   when an exact-match coding DNA sequence is provided.
2. Mature-CDS motif hygiene: eliminates cryptic splice donors ([AC]AGGT[AG]AGT),
   premature polyA signals (AATAAA, ATTAAA, AATAAT), destabilizing AU-rich
   elements (ATTTA), Type IIS sites, and homopolymers (>4 nt).
3. Balanced Arg (AGA:AGG ~ 1:1) and Ser (TCT/TCC/TCA/AGT/AGC ~ 1:1:1:1:1)
   codon distribution across mutable positions.
4. Adjacent identical Arg/Ser tandem-repeat elimination.

All edits are strictly synonymous and verified against GC bounds (40-48% default)
and forbidden motifs before committing.
"""

from __future__ import annotations

import math
import re
from typing import Iterable

from Bio.Data import CodonTable
from Bio.Seq import Seq

ARG_CODONS: tuple[str, ...] = ("AGA", "AGG")
SER_CODONS: tuple[str, ...] = ("TCT", "TCC", "TCA", "AGT", "AGC")
STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})

SPLICE_DONOR_RE = re.compile(r"[AC]AGGT[AG]AGT")
MATURE_RNA_MOTIFS: tuple[str, ...] = ("AATAAA", "ATTAAA", "AATAAT", "ATTTA")

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_AA_TO_CODONS: dict[str, tuple[str, ...]] = {}
for _codon, _aa in CodonTable.standard_dna_table.forward_table.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)  # type: ignore[arg-type]
_AA_TO_CODONS = {aa: tuple(sorted(codons)) for aa, codons in _AA_TO_CODONS.items()}


def _has_forbidden_or_homopolymer(
    seq: str,
    motifs: tuple[str, ...],
    max_homopolymer: int = 4,
    check_mature_rna_motifs: bool = True,
) -> bool:
    for m in motifs:
        if m in seq:
            return True
        rc = m.translate(_COMPLEMENT)[::-1]
        if rc in seq:
            return True
    if check_mature_rna_motifs:
        for m in MATURE_RNA_MOTIFS:
            if m in seq:
                return True
        if SPLICE_DONOR_RE.search(seq):
            return True
    for base in ("A", "C", "G", "T"):
        if base * (max_homopolymer + 1) in seq:
            return True
    return False


def _count_local_violations(seq: str, motifs: tuple[str, ...], max_homopolymer: int = 4) -> int:
    hits = 0
    for m in motifs:
        hits += seq.count(m)
        rc = m.translate(_COMPLEMENT)[::-1]
        if rc != m:
            hits += seq.count(rc)
    for m in MATURE_RNA_MOTIFS:
        hits += seq.count(m)
    hits += len(SPLICE_DONOR_RE.findall(seq))
    for base in ("A", "C", "G", "T"):
        hits += len(re.findall(f"{base}{{{max_homopolymer + 1},}}", seq))
    return hits


def _gc_within_bounds(gc_count: int, total_nt: int, gc_min: float, gc_max: float) -> bool:
    if total_nt <= 0:
        return True
    lo = math.ceil((gc_min / 100.0) * total_nt)
    hi = math.floor((gc_max / 100.0) * total_nt)
    return lo <= gc_count <= hi


def _build_balanced_targets(count: int, palette: tuple[str, ...]) -> dict[str, int]:
    k = len(palette)
    base, rem = divmod(count, k)
    return {codon: base + (1 if idx < rem else 0) for idx, codon in enumerate(palette)}


def apply_default_sop_refinement(
    cds: str,
    *,
    original_dna: str | None = None,
    preserve_signal_peptide: bool = True,
    signal_peptide_nt_len: int = 60,
    gc_min: float = 40.0,
    gc_max: float = 48.0,
    forbidden_motifs: Iterable[str] = (),
    protected_codon_indices: set[int] | None = None,
) -> str:
    """Refine a valid CDS with 5' SP lock (if DNA) and Arg/Ser balancing + repeat guard."""
    clean_cds = "".join(cds.upper().split())
    if len(clean_cds) < 3 or len(clean_cds) % 3 != 0:
        return clean_cds

    has_terminal_stop = clean_cds[-3:] in STOP_CODONS
    stop_tail = clean_cds[-3:] if has_terminal_stop else ""
    body = clean_cds[:-3] if has_terminal_stop else clean_cds
    if not body:
        return clean_cds

    codons = [body[i : i + 3] for i in range(0, len(body), 3)]
    motifs_tuple = tuple({m.upper() for m in forbidden_motifs if m})
    protected = set(protected_codon_indices or ())
    locked_prefix_codons = 0

    # 1. Preserve 5' signal peptide (first 60 nt) when DNA input matches translation
    if (
        preserve_signal_peptide
        and original_dna
        and signal_peptide_nt_len > 0
        and signal_peptide_nt_len % 3 == 0
    ):
        clean_orig = "".join(original_dna.upper().split())
        sp_codons_count = signal_peptide_nt_len // 3
        if len(clean_orig) >= signal_peptide_nt_len and len(codons) >= sp_codons_count:
            orig_prefix = clean_orig[:signal_peptide_nt_len]
            if set(orig_prefix) <= {"A", "C", "G", "T"}:
                orig_aa = str(Seq(orig_prefix).translate())
                curr_aa = str(Seq("".join(codons[:sp_codons_count])).translate())
                if orig_aa == curr_aa:
                    candidate_codons = [
                        orig_prefix[i : i + 3] for i in range(0, signal_peptide_nt_len, 3)
                    ] + codons[sp_codons_count:]
                    candidate_full = "".join(candidate_codons) + stop_tail
                    mature_junction = candidate_full[max(0, signal_peptide_nt_len - 4) :]
                    if not _has_forbidden_or_homopolymer(
                        candidate_full, motifs_tuple, max_homopolymer=999, check_mature_rna_motifs=False
                    ) and not _has_forbidden_or_homopolymer(mature_junction, ()):
                        codons = candidate_codons
                        locked_prefix_codons = sp_codons_count

    # Helper to safely test a single-codon replacement
    total_nt = len(body) + len(stop_tail)
    current_gc = ("".join(codons) + stop_tail).count("G") + ("".join(codons) + stop_tail).count("C")

    def _can_replace(idx: int, new_codon: str, enforce_gc: bool = True) -> bool:
        old_codon = codons[idx]
        if old_codon == new_codon:
            return True
        if enforce_gc:
            lo = math.ceil((gc_min / 100.0) * total_nt)
            hi = math.floor((gc_max / 100.0) * total_nt)
            new_gc = current_gc - (old_codon.count("G") + old_codon.count("C")) + (
                new_codon.count("G") + new_codon.count("C")
            )
            if lo <= current_gc <= hi:
                if not (lo <= new_gc <= hi):
                    return False
            elif current_gc < lo and new_gc < current_gc:
                return False
            elif current_gc > hi and new_gc > current_gc:
                return False
        # Check local window for motifs / homopolymers
        win_start = max(0, idx - 4)
        win_end = min(len(codons), idx + 5)
        snippet = "".join(codons[win_start:idx]) + new_codon + "".join(codons[idx + 1 : win_end])
        if win_end == len(codons) and stop_tail:
            snippet += stop_tail
        return not _has_forbidden_or_homopolymer(snippet, motifs_tuple)

    def _commit(idx: int, new_codon: str) -> None:
        nonlocal current_gc
        old_codon = codons[idx]
        current_gc = current_gc - (old_codon.count("G") + old_codon.count("C")) + (
            new_codon.count("G") + new_codon.count("C")
        )
        codons[idx] = new_codon

    # 2. Eliminate pre-existing mature-CDS SOP motif / homopolymer violations
    body_aa = str(Seq("".join(codons)).translate())
    for i, aa in enumerate(body_aa):
        if i < locked_prefix_codons or i in protected:
            continue
        win_start = max(locked_prefix_codons, i - 4)
        win_end = min(len(codons), i + 5)
        cur_snippet = "".join(codons[win_start:win_end]) + (
            stop_tail if win_end == len(codons) and stop_tail else ""
        )
        cur_hits = _count_local_violations(cur_snippet, motifs_tuple)
        if cur_hits == 0:
            continue
        best_cand = None
        best_hits = cur_hits
        lo_gc = math.ceil((gc_min / 100.0) * total_nt)
        hi_gc = math.floor((gc_max / 100.0) * total_nt)
        for cand in _AA_TO_CODONS.get(aa, ()):
            if cand == codons[i]:
                continue
            new_gc = current_gc - (codons[i].count("G") + codons[i].count("C")) + (
                cand.count("G") + cand.count("C")
            )
            if lo_gc <= current_gc <= hi_gc and not (lo_gc <= new_gc <= hi_gc):
                continue
            if current_gc < lo_gc and new_gc < current_gc:
                continue
            if current_gc > hi_gc and new_gc > current_gc:
                continue
            new_snippet = (
                "".join(codons[win_start:i])
                + cand
                + "".join(codons[i + 1 : win_end])
                + (stop_tail if win_end == len(codons) and stop_tail else "")
            )
            new_hits = _count_local_violations(new_snippet, motifs_tuple)
            if new_hits < best_hits:
                best_hits = new_hits
                best_cand = cand
                if best_hits == 0:
                    break
        if best_cand is not None:
            _commit(i, best_cand)

    # 3. Balance Arg (AGA/AGG) and Ser (5 codons) + avoid adjacent identical repeats
    aa_map = {"R": ARG_CODONS, "S": SER_CODONS}

    for target_aa, palette in aa_map.items():
        mutable_indices = [
            i
            for i, aa in enumerate(body_aa)
            if aa == target_aa and i >= locked_prefix_codons and i not in protected
        ]
        if not mutable_indices:
            continue

        targets = _build_balanced_targets(len(mutable_indices), palette)
        counts = {c: 0 for c in palette}
        for i in mutable_indices:
            if codons[i] in counts:
                counts[codons[i]] += 1

        # Pass 1: Balance excess counts without creating adjacent repeats
        for i in mutable_indices:
            cur = codons[i]
            if cur in targets and counts.get(cur, 0) <= targets[cur]:
                continue

            candidates = sorted(
                palette,
                key=lambda c: (
                    (i > 0 and codons[i - 1] == c) or (i + 1 < len(codons) and codons[i + 1] == c),
                    -(targets[c] - counts[c]),
                    c,
                ),
            )
            for cand in candidates:
                if cand == cur or counts[cand] >= targets[cand]:
                    continue
                if (i > 0 and codons[i - 1] == cand) or (i + 1 < len(codons) and codons[i + 1] == cand):
                    continue
                if _can_replace(i, cand, enforce_gc=True):
                    if cur in counts:
                        counts[cur] -= 1
                    counts[cand] += 1
                    _commit(i, cand)
                    break

        # Pass 2: Eliminate any remaining adjacent identical Arg/Ser tandem repeats
        for i in mutable_indices:
            cur = codons[i]
            adjacent_repeat = (i > 0 and codons[i - 1] == cur) or (
                i + 1 < len(codons) and codons[i + 1] == cur
            )
            if not adjacent_repeat:
                continue
            candidates = sorted(
                palette,
                key=lambda c: (
                    (i > 0 and codons[i - 1] == c) or (i + 1 < len(codons) and codons[i + 1] == c),
                    -(targets[c] - counts[c]),
                    c,
                ),
            )
            for cand in candidates:
                if cand == cur:
                    continue
                if (i > 0 and codons[i - 1] == cand) or (i + 1 < len(codons) and codons[i + 1] == cand):
                    continue
                if _can_replace(i, cand, enforce_gc=True):
                    if cur in counts:
                        counts[cur] -= 1
                    counts[cand] += 1
                    _commit(i, cand)
                    break

    result = "".join(codons) + stop_tail
    # Final safety verification
    if len(result) != len(clean_cds) or str(Seq(result).translate()) != str(Seq(clean_cds).translate()):
        return clean_cds
    return result
