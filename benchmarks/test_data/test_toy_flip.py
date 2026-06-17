"""Toy codon-weight injection acceptance check for benchmark greedy CAI."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
for path in (ROOT, SRC):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from benchmarks.baselines.greedy_cai import greedy_cai_cds
from factorforge.analysis.metrics import load_codon_usage_table


PROTEIN = "MKKKKKFFFY"
PROFILE_A = Path(__file__).with_name("toy_profile_A.json")
PROFILE_B = Path(__file__).with_name("toy_profile_B.json")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _lys_codon_counts(protein: str, cds: str) -> tuple[int, int]:
    codons = [cds[index : index + 3] for index in range(0, len(cds), 3)]
    lys_codons = [codon for aa, codon in zip(protein, codons) if aa == "K"]
    return lys_codons.count("AAA"), lys_codons.count("AAG")


def _load_weights(path: Path) -> dict[str, float]:
    return load_codon_usage_table(path).codon_weights


def main() -> int:
    weights_a = _load_weights(PROFILE_A)
    weights_b = _load_weights(PROFILE_B)

    output_a = greedy_cai_cds(PROTEIN, codon_weights=weights_a)
    output_b = greedy_cai_cds(PROTEIN, codon_weights=weights_b)

    aaa_a, aag_a = _lys_codon_counts(PROTEIN, output_a)
    aaa_b, aag_b = _lys_codon_counts(PROTEIN, output_b)
    sha_a = _sha256(PROFILE_A)
    sha_b = _sha256(PROFILE_B)

    passed = aaa_a > aag_a and aag_b > aaa_b and output_a != output_b
    status = "PASS" if passed else "FAIL"

    print("P2 Toy-Table Flip Acceptance Test")
    print(f"status: {status}")
    print()
    print("| profile | AAA_weight | AAG_weight | output | AAA_count | AAG_count | sha256 |")
    print("|---|---:|---:|---|---:|---:|---|")
    print(
        f"| Toy A | {weights_a['AAA']:.4f} | {weights_a['AAG']:.4f} | "
        f"{output_a} | {aaa_a} | {aag_a} | {sha_a} |"
    )
    print(
        f"| Toy B | {weights_b['AAA']:.4f} | {weights_b['AAG']:.4f} | "
        f"{output_b} | {aaa_b} | {aag_b} | {sha_b} |"
    )
    print()
    print(f"output_A != output_B: {output_a != output_b}")
    print(f"Toy A Lys dominance: AAA_count > AAG_count -> {aaa_a > aag_a}")
    print(f"Toy B Lys dominance: AAG_count > AAA_count -> {aag_b > aaa_b}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
