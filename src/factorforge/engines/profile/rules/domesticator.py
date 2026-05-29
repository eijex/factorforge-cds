"""
Domesticator for FactorForge profile engine
Assembly standard compatibility (P0-4) - Golden Gate/MoClo/BioBricks
"""

from __future__ import annotations

from typing import Any

from factorforge.engines.profile.utils import build_aa_to_codons_map, get_data_path, load_codon_table


class Domesticator:
    """
    Remove restriction enzyme sites for assembly compatibility

    Supported assembly systems:
    - Golden Gate (BsaI, BpiI, BsmBI)
    - MoClo (BsaI + overhangs)
    - BioBricks (EcoRI, XbaI, SpeI, PstI)
    """

    # Assembly standard definitions
    ASSEMBLY_STANDARDS: dict[str, dict[str, Any]] = {
        "golden_gate": {
            "enzymes": ["BsaI", "BpiI", "BsmBI"],
            "sites": {
                "BsaI": ["GGTCTC", "GAGACC"],  # Forward and reverse complement
                "BpiI": ["GAAGAC", "GTCTTC"],
                "BsmBI": ["CGTCTC", "GAGACG"],
            },
        },
        "moclo": {
            "enzymes": ["BsaI"],
            "sites": {"BsaI": ["GGTCTC", "GAGACC"]},
            "overhangs": ["AATG", "AGGT", "GCTT", "CGCT"],  # Level 0
        },
        "biobricks": {
            "enzymes": ["EcoRI", "XbaI", "SpeI", "PstI"],
            "sites": {
                "EcoRI": ["GAATTC"],
                "XbaI": ["TCTAGA"],
                "SpeI": ["ACTAGT"],
                "PstI": ["CTGCAG"],
            },
        },
    }

    def __init__(
        self,
        codon_table: dict[str, Any] | None = None,
        host: str = "nbenthamiana",
    ) -> None:
        """
        Args:
            codon_table: Codon table (loads default if None)
            host: Host codon table name used when codon_table is not provided.
        """
        self.host = host
        if codon_table is None:
            codon_table = load_codon_table(host, get_data_path())

        self.codon_table: dict[str, Any] = codon_table
        self.aa_to_codons: dict[str, list[str]] = self._build_aa_to_codons_map()

    def _build_aa_to_codons_map(self) -> dict[str, list[str]]:
        """Build amino-acid-to-codons map"""
        return build_aa_to_codons_map(self.codon_table)

    def scan_restriction_sites(
        self,
        seq: str,
        standard: str = "golden_gate",
    ) -> list[dict[str, Any]]:
        """
        Scan restriction enzyme sites

        Args:
            seq: DNA sequence
            standard: Assembly standard ("golden_gate", "moclo", "biobricks")

        Returns:
            List of detected sites

        Raises:
            ValueError: Unsupported assembly standard.

        Examples:
            >>> domesticator = Domesticator()
            >>> domesticator.scan_restriction_sites("GGTCTC", "golden_gate")
            [{'enzyme': 'BsaI', ...}]
        """
        if standard not in self.ASSEMBLY_STANDARDS:
            raise ValueError(f"Unknown assembly standard: {standard}")

        assembly_info = self.ASSEMBLY_STANDARDS[standard]
        sites_found: list[dict[str, Any]] = []

        for enzyme, site_seqs in assembly_info["sites"].items():
            for site_seq in site_seqs:
                pos = 0
                while True:
                    idx = seq.find(site_seq, pos)
                    if idx == -1:
                        break

                    sites_found.append(
                        {
                            "enzyme": enzyme,
                            "site": site_seq,
                            "position": idx,
                            "context": seq[
                                max(0, idx - 10) : min(len(seq), idx + len(site_seq) + 10)
                            ],
                        }
                    )
                    pos = idx + 1

        return sites_found

    def domesticate(
        self,
        seq: str,
        standard: str = "golden_gate",
        max_attempts: int = 100,
    ) -> dict[str, Any]:
        """
        Remove restriction enzyme sites

        Args:
            seq: DNA sequence
            standard: Assembly standard
            max_attempts: Maximum attempts

        Returns:
            {
                "domesticated_seq": "...",
                "removed_sites": [...],
                "unfixable": [...],
                "success": True/False
            }

        Raises:
            ValueError: Unsupported assembly standard.

        Examples:
            >>> domesticator = Domesticator()
            >>> result = domesticator.domesticate("ATGGGTCTCGAG", "golden_gate")
            >>> "domesticated_seq" in result
            True
        """
        if len(seq) % 3 != 0:
            return {
                "success": False,
                "error": "Sequence length not divisible by 3",
                "domesticated_seq": seq,
            }

        modified_seq = seq
        removed_sites: list[dict[str, Any]] = []
        unfixable: list[dict[str, Any]] = []

        # Iteratively remove sites
        for attempt in range(max_attempts):
            sites = self.scan_restriction_sites(modified_seq, standard)

            if not sites:
                # All sites removed
                return {
                    "success": True,
                    "domesticated_seq": modified_seq,
                    "removed_sites": removed_sites,
                    "unfixable": [],
                    "attempts": attempt + 1,
                }

            # Attempt to remove the first site
            site = sites[0]
            result = self._remove_site(modified_seq, site)

            if result["success"]:
                modified_seq = result["modified_seq"]
                removed_sites.append(
                    {
                        "enzyme": site["enzyme"],
                        "position": site["position"],
                        "changes": result["changes"],
                    }
                )
            else:
                # Failed to remove site
                unfixable.append(
                    {
                        "enzyme": site["enzyme"],
                        "site": site["site"],
                        "position": site["position"],
                        "reason": result.get("reason", "Unknown"),
                        "alternatives": self._suggest_alternatives(seq, site),
                    }
                )
                # Record as unfixable and continue to next site
                continue

        # Sites still remain after max attempts
        remaining_sites = self.scan_restriction_sites(modified_seq, standard)

        return {
            "success": len(remaining_sites) == 0,
            "domesticated_seq": modified_seq,
            "removed_sites": removed_sites,
            "unfixable": unfixable if unfixable else remaining_sites,
            "attempts": max_attempts,
        }

    def _remove_site(self, seq: str, site: dict[str, Any]) -> dict[str, Any]:
        """
        Remove a single restriction enzyme site

        Args:
            seq: DNA sequence
            site: Site info

        Returns:
            {
                "success": True/False,
                "modified_seq": "...",
                "changes": [...]
            }
        """
        pos = site["position"]
        site_seq = site["site"]
        site_len = len(site_seq)

        # Compute codon range overlapping the site
        first_codon_idx = (pos // 3) * 3
        last_codon_idx = ((pos + site_len - 1) // 3) * 3

        # Try synonymous substitutions per codon
        for codon_start in range(first_codon_idx, last_codon_idx + 1, 3):
            if codon_start + 3 > len(seq):
                continue

            original_codon = seq[codon_start : codon_start + 3]

            # Validate amino acid
            if original_codon not in self.codon_table["codons"]:
                continue

            aa = self.codon_table["codons"][original_codon]["aa"]

            # Find synonymous codons
            synonymous_codons = [c for c in self.aa_to_codons.get(aa, []) if c != original_codon]

            if not synonymous_codons:
                continue

            # Try each synonymous codon
            for alt_codon in synonymous_codons:
                # Temporary substitution
                test_seq = seq[:codon_start] + alt_codon + seq[codon_start + 3 :]

                # Check if site is gone
                test_region = test_seq[max(0, pos - 10) : min(len(test_seq), pos + site_len + 10)]

                if site_seq not in test_region:
                    # Success
                    return {
                        "success": True,
                        "modified_seq": test_seq,
                        "changes": [
                            {
                                "pos": codon_start,
                                "original": original_codon,
                                "fixed": alt_codon,
                                "aa": aa,
                            }
                        ],
                    }

        # Failed to fix
        return {
            "success": False,
            "modified_seq": seq,
            "changes": [],
            "reason": "No synonymous codon available to remove site",
        }

    def _suggest_alternatives(self, seq: str, site: dict[str, Any]) -> list[str]:
        """
        Suggest alternatives for unremovable sites

        Args:
            seq: DNA sequence
            site: Site info

        Returns:
            List of alternative suggestions
        """
        alternatives: list[str] = []

        pos = site["position"]
        site_len = len(site["site"])

        # Determine affected codons
        first_codon_idx = (pos // 3) * 3
        last_codon_idx = ((pos + site_len - 1) // 3) * 3

        affected_codons: list[dict[str, Any]] = []
        for codon_start in range(first_codon_idx, last_codon_idx + 1, 3):
            if codon_start + 3 <= len(seq):
                codon = seq[codon_start : codon_start + 3]
                if codon in self.codon_table["codons"]:
                    aa = self.codon_table["codons"][codon]["aa"]
                    synonymous = self.aa_to_codons.get(aa, [])
                    affected_codons.append(
                        {
                            "pos": codon_start,
                            "codon": codon,
                            "aa": aa,
                            "synonymous_count": len(synonymous) - 1,
                        }
                    )

        # Suggest alternatives
        if any(c["synonymous_count"] == 0 for c in affected_codons):
            alternatives.append(
                "Includes amino acids without synonyms - requires non-synonymous change"
            )

        alternatives.append("Try shifting the site by adjusting adjacent codons")

        alternatives.append("Consider using a different assembly method")

        return alternatives

    def batch_domesticate(
        self,
        sequences: list[dict[str, Any]],
        standard: str = "golden_gate",
    ) -> list[dict[str, Any]]:
        """
        Batch domestication

        Args:
            sequences: [{"id": "gene1", "sequence": "ATG..."}, ...]
            standard: Assembly standard

        Returns:
            List of results

        Raises:
            ValueError: Unsupported assembly standard.

        Examples:
            >>> domesticator = Domesticator()
            >>> results = domesticator.batch_domesticate([{"id": "x", "sequence": "ATG"}])
            >>> len(results)
            1
        """
        results: list[dict[str, Any]] = []

        for seq_data in sequences:
            seq_id = seq_data.get("id", "unknown")
            seq = seq_data.get("sequence", "")

            result = self.domesticate(seq, standard)
            result["id"] = seq_id
            results.append(result)

        return results


# --- Usage example ---
if __name__ == "__main__":
    domesticator = Domesticator()

    # Test sequence (includes BsaI site)
    test_seq = "ATGGGTCTCGAGGAGCTGTTCACCGGGGTGGTGCCCATC"

    print("=== Original Sequence ===")
    print(f"Sequence: {test_seq}")

    # Golden Gate scan
    sites = domesticator.scan_restriction_sites(test_seq, "golden_gate")
    print(f"\nRestriction sites found: {len(sites)}")
    for site in sites:
        print(f"  - {site['enzyme']} at position {site['position']}: {site['site']}")

    # Domestication
    print("\n=== Domestication ===")
    result = domesticator.domesticate(test_seq, "golden_gate")

    if result["success"]:
        print("✅ Domestication successful!")
        print(f"Modified sequence: {result['domesticated_seq']}")
        print(f"Removed sites: {len(result['removed_sites'])}")
        for removed in result["removed_sites"]:
            print(f"  - {removed['enzyme']} at position {removed['position']}")
    else:
        print("❌ Domestication failed")
        print(f"Unfixable sites: {len(result['unfixable'])}")
        for unfixable in result["unfixable"]:
            print(
                f"  - {unfixable.get('enzyme', 'N/A')} at position {unfixable.get('position', 'N/A')}"
            )
            print(f"    Alternatives: {unfixable.get('alternatives', [])}")
