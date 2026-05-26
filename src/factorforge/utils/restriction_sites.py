"""Custom restriction-site detection and CDS domestication utilities."""

from __future__ import annotations

from typing import Any, Iterable

from factorforge.engines.profile.utils import build_aa_to_codons_map

DNA_BASES = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}
_RC_TRANS = str.maketrans("ACGT", "TGCA")


BUILT_IN_RESTRICTION_SITES: tuple[dict[str, Any], ...] = (
    {"name": "BsaI", "sequence": "GGTCTC", "scan_rc": True},
    {"name": "BpiI", "sequence": "GAAGAC", "scan_rc": True},
    {"name": "BsmBI", "sequence": "CGTCTC", "scan_rc": True},
)


def detect_restriction_sites(
    seq: str,
    sites: Iterable[dict[str, Any] | tuple[str, str] | str],
    scan_rc: bool = True,
) -> list[dict[str, Any]]:
    """Detect custom restriction sites in a DNA sequence.

    Args:
        seq: DNA sequence to scan. Whitespace is ignored.
        sites: Site definitions as dicts with ``name`` and ``sequence``, as
            ``(name, sequence)`` tuples, or as sequence strings.
        scan_rc: Default reverse-complement scanning behavior.

    Returns:
        List of hits with ``name``, ``sequence``, ``position``, and ``strand``.

    Raises:
        ValueError: If the sequence or site definitions are invalid.
    """
    normalized_seq = _normalize_dna(seq, field_name="sequence", allow_empty=True)
    normalized_sites = _normalize_sites(sites, default_scan_rc=scan_rc)
    return _public_hits(_detect_with_normalized_sites(normalized_seq, normalized_sites))


def remove_restriction_site(
    seq: str,
    hit: dict[str, Any],
    codon_table: dict[str, Any],
) -> tuple[str, str] | None:
    """Remove one restriction-site hit by synonymous codon substitution.

    The function preserves amino-acid identity and refuses substitutions that
    introduce any built-in Golden Gate site or any site supplied through the
    private ``_all_sites`` hit key used by ``domesticate_custom_sites``.
    """
    normalized_seq = _normalize_dna(seq, field_name="sequence", allow_empty=True)
    codons_section = _codons_section(codon_table)
    pos = int(hit["position"])
    site_sequence = _normalize_dna(hit["sequence"], field_name="site sequence")
    matched_sequence = (
        _reverse_complement(site_sequence) if hit.get("strand") == "rc" else site_sequence
    )
    site_end = pos + len(matched_sequence)
    editable_end = _editable_cds_end(normalized_seq)

    if pos < 0 or site_end > editable_end:
        return None

    aa_to_codons = build_aa_to_codons_map(codon_table)
    before_translation = _translate_with_table(normalized_seq, codon_table).rstrip("*")
    all_sites = hit.get("_all_sites") or (
        list(BUILT_IN_RESTRICTION_SITES)
        + [{"name": hit["name"], "sequence": site_sequence, "scan_rc": True}]
    )
    before_site_keys = _site_keys(_detect_with_normalized_sites(normalized_seq, all_sites))

    first_codon_idx = (pos // 3) * 3
    last_codon_idx = ((site_end - 1) // 3) * 3

    for codon_start in range(first_codon_idx, last_codon_idx + 1, 3):
        if codon_start + 3 > editable_end:
            continue

        original_codon = normalized_seq[codon_start : codon_start + 3]
        codon_info = codons_section.get(original_codon)
        if not codon_info:
            continue

        aa = str(codon_info.get("aa", ""))
        if aa == "*":
            continue

        alternatives = _synonymous_codons(
            aa=aa,
            original_codon=original_codon,
            aa_to_codons=aa_to_codons,
            codons_section=codons_section,
        )
        for alt_codon in alternatives:
            candidate = normalized_seq[:codon_start] + alt_codon + normalized_seq[codon_start + 3 :]
            if candidate[pos:site_end] == matched_sequence:
                continue
            if not _strict_cds_guard(candidate, codon_table, before_translation):
                continue

            after_site_keys = _site_keys(_detect_with_normalized_sites(candidate, all_sites))
            if after_site_keys - before_site_keys:
                continue

            return candidate, f"{original_codon}->{alt_codon}"

    return None


def domesticate_custom_sites(
    seq: str,
    sites: Iterable[dict[str, Any] | tuple[str, str] | str],
    codon_table: dict[str, Any],
    max_iter: int = 10,
) -> dict[str, Any]:
    """Iteratively remove custom restriction sites from a CDS.

    Sites outside the editable CDS region are reported in ``unresolved`` and
    left unchanged. If the input CDS has a terminal stop codon, that final codon
    is treated as non-editable.
    """
    if max_iter < 1:
        raise ValueError("max_iter must be >= 1")

    modified_seq = _normalize_dna(seq, field_name="sequence", allow_empty=True)
    custom_sites = _normalize_sites(sites, default_scan_rc=True)
    all_sites = [*custom_sites, *BUILT_IN_RESTRICTION_SITES]
    detected = _detect_with_normalized_sites(modified_seq, custom_sites)
    removed: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    unresolved_keys: set[tuple[str, int, str]] = set()
    iterations = 0

    if not _strict_cds_guard(
        modified_seq,
        codon_table,
        _translate_with_table(modified_seq, codon_table).rstrip("*"),
    ):
        raise ValueError("Input sequence failed strict CDS validation")

    for iteration in range(1, max_iter + 1):
        iterations = iteration
        current_hits = _detect_with_normalized_sites(modified_seq, custom_sites)
        pending_hits = [hit for hit in current_hits if _unresolved_key(hit) not in unresolved_keys]

        if not pending_hits:
            break

        hit = pending_hits[0]
        hit["_all_sites"] = all_sites

        if _is_outside_editable_cds(modified_seq, hit):
            unresolved.append(_unresolved_hit(hit, reason="outside_cds"))
            unresolved_keys.add(_unresolved_key(hit))
            continue

        removal = remove_restriction_site(modified_seq, hit, codon_table)
        if removal is None:
            unresolved.append(_unresolved_hit(hit, reason="no_synonymous_substitution"))
            unresolved_keys.add(_unresolved_key(hit))
            continue

        modified_seq, substitution = removal
        removed.append(
            {
                "name": hit["name"],
                "position": hit["position"],
                "substitution": substitution,
            }
        )

    remaining_hits = _detect_with_normalized_sites(modified_seq, custom_sites)
    for hit in remaining_hits:
        key = _unresolved_key(hit)
        if key not in unresolved_keys:
            unresolved.append(_unresolved_hit(hit, reason="max_iter_reached"))
            unresolved_keys.add(key)

    if not _strict_cds_guard(
        modified_seq,
        codon_table,
        _translate_with_table(_normalize_dna(seq, field_name="sequence"), codon_table).rstrip("*"),
    ):
        raise ValueError("Domesticated sequence failed strict CDS validation")

    return {
        "sequence": modified_seq,
        "detected": _public_hits(detected),
        "removed": removed,
        "unresolved": unresolved,
        "iterations": iterations,
    }


def _normalize_dna(seq: str, field_name: str, allow_empty: bool = False) -> str:
    normalized = "".join(str(seq).upper().split()).replace("U", "T")
    if not normalized and not allow_empty:
        raise ValueError(f"{field_name} must not be empty")
    invalid = set(normalized) - DNA_BASES
    if invalid:
        chars = ", ".join(sorted(invalid))
        raise ValueError(f"{field_name} must contain only A/C/G/T bases: {chars}")
    return normalized


def _normalize_sites(
    sites: Iterable[dict[str, Any] | tuple[str, str] | str],
    default_scan_rc: bool,
) -> list[dict[str, Any]]:
    normalized_sites: list[dict[str, Any]] = []
    seen_names: set[str] = set()
    seen_sequences: set[str] = set()

    for raw_site in sites:
        if isinstance(raw_site, dict):
            raw_sequence = raw_site.get("sequence", raw_site.get("site", ""))
            sequence = _normalize_dna(raw_sequence, field_name="site sequence")
            name = str(raw_site.get("name") or sequence).strip()
            site_scan_rc = bool(raw_site.get("scan_rc", default_scan_rc))
        elif isinstance(raw_site, tuple) and len(raw_site) == 2:
            name = str(raw_site[0]).strip()
            sequence = _normalize_dna(raw_site[1], field_name="site sequence")
            site_scan_rc = default_scan_rc
        else:
            sequence = _normalize_dna(str(raw_site), field_name="site sequence")
            name = sequence
            site_scan_rc = default_scan_rc

        if not name:
            raise ValueError("site name must not be empty")
        if len(sequence) < 4 or len(sequence) > 12:
            raise ValueError("site sequence length must be 4-12 bp")
        if name in seen_names:
            raise ValueError(f"duplicate site name: {name}")
        if sequence in seen_sequences:
            raise ValueError(f"duplicate site sequence: {sequence}")

        seen_names.add(name)
        seen_sequences.add(sequence)
        normalized_sites.append({"name": name, "sequence": sequence, "scan_rc": site_scan_rc})

    return normalized_sites


def _detect_with_normalized_sites(
    seq: str,
    sites: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    for site in sites:
        patterns = [("fwd", str(site["sequence"]))]
        rc_sequence = _reverse_complement(str(site["sequence"]))
        if bool(site.get("scan_rc", True)) and rc_sequence != site["sequence"]:
            patterns.append(("rc", rc_sequence))

        for strand, pattern in patterns:
            pos = 0
            while True:
                idx = seq.find(pattern, pos)
                if idx == -1:
                    break
                hits.append(
                    {
                        "name": site["name"],
                        "sequence": site["sequence"],
                        "position": idx,
                        "strand": strand,
                    }
                )
                pos = idx + 1

    return sorted(hits, key=lambda item: (int(item["position"]), str(item["name"]), item["strand"]))


def _public_hits(hits: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "name": hit["name"],
            "sequence": hit["sequence"],
            "position": hit["position"],
            "strand": hit["strand"],
        }
        for hit in hits
    ]


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TRANS)[::-1]


def _editable_cds_end(seq: str) -> int:
    if len(seq) >= 3 and seq[-3:] in STOP_CODONS:
        return len(seq) - 3
    return len(seq)


def _is_outside_editable_cds(seq: str, hit: dict[str, Any]) -> bool:
    pos = int(hit["position"])
    site_len = len(str(hit["sequence"]))
    return pos < 0 or pos + site_len > _editable_cds_end(seq)


def _codons_section(codon_table: dict[str, Any]) -> dict[str, dict[str, Any]]:
    codons = codon_table.get("codons")
    if not isinstance(codons, dict):
        raise ValueError("codon_table must contain a codons mapping")
    return codons


def _synonymous_codons(
    aa: str,
    original_codon: str,
    aa_to_codons: dict[str, list[str]],
    codons_section: dict[str, dict[str, Any]],
) -> list[str]:
    codons = [codon for codon in aa_to_codons.get(aa, []) if codon != original_codon]
    return sorted(
        codons,
        key=lambda codon: float(codons_section.get(codon, {}).get("frequency", 0.0)),
        reverse=True,
    )


def _translate_with_table(seq: str, codon_table: dict[str, Any]) -> str:
    codons_section = _codons_section(codon_table)
    translated: list[str] = []
    for index in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[index : index + 3]
        translated.append(str(codons_section.get(codon, {}).get("aa", "X")))
    return "".join(translated)


def _strict_cds_guard(seq: str, codon_table: dict[str, Any], expected_translation: str) -> bool:
    if len(seq) % 3 != 0:
        return False

    codons_section = _codons_section(codon_table)
    codons = [seq[index : index + 3] for index in range(0, len(seq), 3)]
    if any(codon not in codons_section for codon in codons):
        return False

    translated = _translate_with_table(seq, codon_table)
    if "*" in translated[:-1]:
        return False

    return translated.rstrip("*") == expected_translation


def _site_keys(hits: Iterable[dict[str, Any]]) -> set[tuple[str, str, int, str]]:
    return {
        (str(hit["name"]), str(hit["sequence"]), int(hit["position"]), str(hit["strand"]))
        for hit in hits
    }


def _unresolved_key(hit: dict[str, Any]) -> tuple[str, int, str]:
    return (str(hit["name"]), int(hit["position"]), str(hit["strand"]))


def _unresolved_hit(hit: dict[str, Any], reason: str) -> dict[str, Any]:
    return {
        "name": hit["name"],
        "position": hit["position"],
        "reason": reason,
    }
