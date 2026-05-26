"""
Construct Builder for FactorForge profile engine.
Builds Golden Gate-compatible expression constructs from templates.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


class ConstructBuilder:
    """Assemble constructs from JSON templates."""

    def __init__(self, template_dir: Path) -> None:
        """
        Args:
            template_dir: Directory containing construct templates.
        """
        self.template_dir = template_dir

    def load_template(self, name: str) -> dict[str, Any]:
        """
        Load a construct template by name.

        Args:
            name: Template name or filename (without extension).

        Returns:
            Template payload as a dictionary.

        Raises:
            FileNotFoundError: If the template file does not exist.
            json.JSONDecodeError: If the template file is invalid JSON.
            ValueError: If the template path is invalid (path traversal attempt).
        """
        filename = name if name.endswith(".json") else f"{name}.json"
        template_path = self.template_dir / filename

        # Security: Prevent path traversal attacks
        try:
            resolved_path = template_path.resolve()
            resolved_dir = self.template_dir.resolve()
            if not resolved_path.is_relative_to(resolved_dir):
                raise ValueError(
                    f"Invalid template path: {name}. "
                    "Template must be within the template directory."
                )
        except (ValueError, OSError) as exc:
            raise ValueError(f"Invalid template path: {name}") from exc

        with open(template_path, "r", encoding="utf-8") as handle:
            return cast(dict[str, Any], json.load(handle))

    def assemble_parts(self, gene_sequence: str, template: dict[str, Any]) -> str:
        """
        Assemble a construct sequence from template components.

        Args:
            gene_sequence: Optimized CDS sequence.
            template: Template dictionary from load_template().

        Returns:
            Assembled construct DNA sequence.

        Raises:
            ValueError: If template components are missing.
        """
        components = template.get("components", [])
        if not components:
            raise ValueError("Template has no components to assemble.")

        parts: list[str] = []
        for component in components:
            sequence = component.get("sequence", "")
            if component.get("type") == "cds" and sequence == "USER_INPUT":
                sequence = gene_sequence
            parts.append(sequence)

        return "".join(parts)

    def add_features(self, construct_seq: str, template: dict[str, Any]) -> "SeqRecord":
        """
        Create a SeqRecord with component features.

        Args:
            construct_seq: Assembled construct sequence.
            template: Template dictionary.

        Returns:
            SeqRecord with component features added.

        Raises:
            ImportError: If Biopython is not installed.
            ValueError: If component lengths cannot be resolved.
        """
        try:
            from Bio.Seq import Seq
            from Bio.SeqFeature import FeatureLocation, SeqFeature
            from Bio.SeqRecord import SeqRecord
        except ImportError as exc:
            raise ImportError("Biopython is required: pip install biopython") from exc

        components: list[dict[str, Any]] = template.get("components", [])
        template_name = template.get("name", "Construct")
        template_desc = template.get("description", "")
        record_id = template_name.replace(" ", "_")
        record_name = record_id[:16]

        record = SeqRecord(
            Seq(construct_seq),
            id=record_id,
            name=record_name,
            description=template_desc,
        )
        record.annotations["molecule_type"] = "DNA"

        lengths: list[int | None] = []
        unknown_indices: list[int] = []
        for idx, component in enumerate(components):
            sequence = component.get("sequence", "")
            if component.get("type") == "cds" and sequence == "USER_INPUT":
                lengths.append(None)
                unknown_indices.append(idx)
            else:
                lengths.append(len(sequence))

        if len(unknown_indices) > 1:
            raise ValueError("Multiple USER_INPUT components are not supported.")

        if unknown_indices:
            known_total = sum(length for length in lengths if length is not None)
            unknown_length = len(construct_seq) - known_total
            if unknown_length < 0:
                raise ValueError("Construct sequence shorter than template components.")
            lengths[unknown_indices[0]] = unknown_length

        feature_type_map = {
            "promoter": "promoter",
            "5utr": "5'UTR",
            "cds": "CDS",
            "terminator": "terminator",
        }

        cursor = 0
        for component, length in zip(components, lengths):
            if length is None:
                raise ValueError("Component length could not be resolved.")
            start = cursor
            end = cursor + length
            cursor = end

            comp_type = component.get("type", "misc_feature")
            feature_type = feature_type_map.get(comp_type, comp_type)
            label = component.get("name", comp_type)

            feature = SeqFeature(  # type: ignore[no-untyped-call]
                FeatureLocation(start, end),  # type: ignore[no-untyped-call]
                type=feature_type,
                qualifiers={
                    "label": [label],
                    "note": [comp_type],
                },
            )
            record.features.append(feature)

        return record

    def validate_construct(
        self, construct: "SeqRecord", template: dict[str, Any]
    ) -> tuple[bool, list[str]]:
        """
        Validate an assembled construct.

        Args:
            construct: SeqRecord with assembled sequence.
            template: Template dictionary.

        Returns:
            Tuple of (valid, warnings).
        """
        warnings: list[str] = []
        valid = True

        seq_str = str(construct.seq)
        seq_len = len(seq_str)

        if seq_len < 500 or seq_len > 20000:
            warnings.append(f"Construct length {seq_len} bp is outside expected range (500-20000).")
            valid = False

        expected_features = len(template.get("components", []))
        actual_features = len(construct.features)
        if actual_features != expected_features:
            warnings.append(
                f"Feature count {actual_features} does not match template ({expected_features})."
            )
            valid = False

        restriction_sites = {
            "BsaI": ["GGTCTC", "GAGACC"],
            "BpiI": ["GAAGAC", "GTCTTC"],
            "BsmBI": ["CGTCTC", "GAGACG"],
        }
        for enzyme, motifs in restriction_sites.items():
            for motif in motifs:
                if motif in seq_str:
                    warnings.append(f"{enzyme} site detected: {motif}")
                    break

        polya_patterns = ["AATAAA", "ATTAAA", "AGTAAA"]
        for feature in construct.features:
            if feature.type != "CDS":
                continue
            if feature.location is None:
                warnings.append("CDS feature has no location defined.")
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            cds_seq = seq_str[start:end]
            for pattern in polya_patterns:
                if pattern in cds_seq:
                    warnings.append(f"PolyA signal {pattern} detected in CDS.")
                    break

            # Check internal overhang collisions within CDS
            collisions = self.check_internal_overhang_collisions(cds_seq)
            for collision in collisions:
                warnings.append(
                    f"MoClo overhang '{collision['overhang']}' found internally in CDS "
                    f"at position {collision['position']} ({collision['strand']})."
                )

        # Positive PolyA check: terminator/3'UTR must contain a PolyA signal
        for feature in construct.features:
            if feature.type != "terminator":
                continue
            if feature.location is None:
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            term_seq = seq_str[start:end]
            has_polya = any(pattern in term_seq for pattern in polya_patterns)
            if not has_polya:
                warnings.append(
                    "No PolyA signal found in terminator region. "
                    "This may impair mRNA polyadenylation."
                )

        return valid, warnings

    def generate_construct(self, gene_sequence: str, template_name: str) -> "SeqRecord":
        """
        Generate a construct from a template name and gene sequence.

        Args:
            gene_sequence: Optimized CDS sequence.
            template_name: Template name (e.g., "standard_expression").

        Returns:
            SeqRecord with features.
        """
        template = self.load_template(template_name)
        construct_seq = self.assemble_parts(gene_sequence, template)
        construct = self.add_features(construct_seq, template)
        valid, warnings = self.validate_construct(construct, template)

        if warnings:
            status = "VALID" if valid else "INVALID"
            log_func = logger.warning if not valid else logger.info
            log_func(f"Construct {status}: {len(warnings)} warning(s)")
            for warning in warnings:
                log_func(f" - {warning}")

        return construct

    # MoClo Level 0 standard overhangs for CDS parts
    MOCLO_LEVEL0_OVERHANGS: dict[str, str] = {
        "cds_5prime": "AATG",
        "cds_3prime": "GCTT",
    }

    def validate_overhangs(
        self,
        parts: list[dict[str, Any]],
        standard: str = "moclo_level0",
    ) -> tuple[bool, list[str]]:
        """
        Validate Golden Gate overhang consistency for ordered parts.

        For MoClo Level 0 CDS standard:
        - 5' overhang must be AATG
        - 3' overhang must be GCTT
        - Adjacent parts: 3' overhang of part N must match 5' overhang of part N+1

        Args:
            parts: Ordered list of part dictionaries with 'overhang_5' and 'overhang_3' keys.
            standard: Assembly standard to validate against.

        Returns:
            Tuple of (valid, warnings).
        """
        warnings: list[str] = []

        if not parts:
            warnings.append("No parts provided for overhang validation.")
            return False, warnings

        if standard == "moclo_level0":
            expected_5 = self.MOCLO_LEVEL0_OVERHANGS["cds_5prime"]
            expected_3 = self.MOCLO_LEVEL0_OVERHANGS["cds_3prime"]

            # Check first part 5' overhang
            first_oh5 = parts[0].get("overhang_5", "")
            if first_oh5 and first_oh5 != expected_5:
                warnings.append(
                    f"First part 5' overhang '{first_oh5}' does not match "
                    f"MoClo Level 0 expected '{expected_5}'."
                )

            # Check last part 3' overhang
            last_oh3 = parts[-1].get("overhang_3", "")
            if last_oh3 and last_oh3 != expected_3:
                warnings.append(
                    f"Last part 3' overhang '{last_oh3}' does not match "
                    f"MoClo Level 0 expected '{expected_3}'."
                )

        # Check chain consistency: part N 3' overhang == part N+1 5' overhang
        for i in range(len(parts) - 1):
            oh3 = parts[i].get("overhang_3", "")
            oh5_next = parts[i + 1].get("overhang_5", "")
            if oh3 and oh5_next and oh3 != oh5_next:
                warnings.append(
                    f"Overhang mismatch between part {i} (3'={oh3}) "
                    f"and part {i + 1} (5'={oh5_next})."
                )

        valid = len(warnings) == 0
        return valid, warnings

    def check_internal_overhang_collisions(
        self,
        cds_seq: str,
        overhangs: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        """
        Check for internal occurrences of MoClo overhang sequences within CDS.

        Scans for both forward and reverse complement of each overhang to prevent
        assembly artifacts during Golden Gate cloning.

        Args:
            cds_seq: Coding DNA sequence to scan.
            overhangs: List of 4bp overhang sequences to check.
                       Defaults to MoClo Level 0 CDS overhangs [AATG, GCTT].

        Returns:
            List of collision dicts with 'overhang', 'position', 'strand' keys.
        """
        if overhangs is None:
            overhangs = list(self.MOCLO_LEVEL0_OVERHANGS.values())

        # Build reverse complement lookup
        complement = str.maketrans("ATGC", "TACG")
        collisions: list[dict[str, Any]] = []

        for overhang in overhangs:
            rc = overhang.translate(complement)[::-1]

            for i in range(len(cds_seq) - len(overhang) + 1):
                fragment = cds_seq[i : i + len(overhang)]
                if fragment == overhang:
                    collisions.append({"overhang": overhang, "position": i, "strand": "forward"})
                elif fragment == rc:
                    collisions.append(
                        {"overhang": overhang, "position": i, "strand": "reverse_complement"}
                    )

        return collisions

    def assemble_construct(self, gene: str, template: dict[str, Any]) -> "SeqRecord":
        """
        Assemble a construct from a gene sequence and template.

        Args:
            gene: Optimized CDS sequence to insert.
            template: Loaded template dictionary.

        Returns:
            SeqRecord for the assembled construct.
        """
        construct_seq = self.assemble_parts(gene, template)
        return self.add_features(construct_seq, template)
