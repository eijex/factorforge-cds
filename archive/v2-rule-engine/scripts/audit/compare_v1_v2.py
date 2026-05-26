from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Sequence, Tuple

from factorforge.engines.v1_archived.inference import InferenceConfig, infer_protein_sequence
from factorforge.engines.v2.rules.domesticator import Domesticator
from factorforge.engines.v2.rules.reverse_translator import OptimizationProfile, ReverseTranslator
from factorforge.engines.v2.rules.rule_engine import RuleEngine
from factorforge.engines.v2.validator import InputValidator


@dataclass(frozen=True)
class ComparisonConfig:
    profile: OptimizationProfile = OptimizationProfile.BALANCED
    assembly_standard: str = "golden_gate"
    window_size_codons: int = 120
    stride_codons: int = 100


def load_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq_parts: List[str] = []

    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip() or "sequence"
            seq_parts = []
        else:
            seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def compare_sequences(
    sequences: Sequence[Tuple[str, str]],
    config: ComparisonConfig,
) -> List[dict[str, Any]]:
    translator = ReverseTranslator()
    rule_engine = RuleEngine()
    domesticator = Domesticator()
    validator = InputValidator()

    results: List[dict[str, Any]] = []
    for seq_id, protein_seq in sequences:
        validation = validator.validate(protein_seq)
        if not validation["valid"]:
            raise ValueError(f"Invalid protein sequence {seq_id}: {validation['errors']}")

        processed_seq = validation["processed_sequence"]

        v1_dna = infer_protein_sequence(
            processed_seq,
            config=InferenceConfig(
                window_size_codons=config.window_size_codons,
                stride_codons=config.stride_codons,
            ),
        )

        v2_candidates = translator.generate_candidates(
            processed_seq, profile=config.profile, n=1
        )
        if not v2_candidates:
            raise ValueError(f"No v2 candidates generated for {seq_id}")
        v2_dna = v2_candidates[0]["sequence"]

        result = {
            "id": seq_id,
            "input_protein": processed_seq,
            "v1": _analyze_dna(
                v1_dna,
                translator=translator,
                rule_engine=rule_engine,
                domesticator=domesticator,
                validator=validator,
                assembly_standard=config.assembly_standard,
            ),
            "v2": _analyze_dna(
                v2_dna,
                translator=translator,
                rule_engine=rule_engine,
                domesticator=domesticator,
                validator=validator,
                assembly_standard=config.assembly_standard,
            ),
        }
        results.append(result)

    return results


def _analyze_dna(
    dna_sequence: str,
    *,
    translator: ReverseTranslator,
    rule_engine: RuleEngine,
    domesticator: Domesticator,
    validator: InputValidator,
    assembly_standard: str,
) -> dict[str, Any]:
    metrics = {
        "length_bp": len(dna_sequence),
        "cai": translator.calculate_cai(dna_sequence),
        "gc": translator.calculate_gc_content(dna_sequence),
        "polya_count": len(rule_engine.scan_polya(dna_sequence)),
        "restriction_sites": len(
            domesticator.scan_restriction_sites(dna_sequence, standard=assembly_standard)
        ),
        "splice_sites": len(rule_engine.scan_splice_sites(dna_sequence)),
    }

    validation = validator.validate(dna_sequence)
    domestication = domesticator.domesticate(dna_sequence, standard=assembly_standard)

    return {
        "sequence": dna_sequence,
        "metrics": metrics,
        "validation": {
            "valid": validation["valid"],
            "level": validation["level"],
            "warnings": len(validation["warnings"]),
            "errors": len(validation["errors"]),
        },
        "domestication": {
            "success": domestication.get("success", False),
            "removed_sites": len(domestication.get("removed_sites", [])),
            "unfixable_sites": len(domestication.get("unfixable", []))
            if isinstance(domestication.get("unfixable"), list)
            else 0,
        },
    }


def render_table(results: Iterable[dict[str, Any]]) -> str:
    lines: List[str] = []
    header = (
        "engine",
        "len_bp",
        "cai",
        "gc",
        "polya",
        "restr_sites",
        "splice",
        "dom_removed",
        "dom_success",
        "valid",
    )
    lines.append(" | ".join(header))
    lines.append("-" * (len(lines[0]) + 2))

    for result in results:
        for engine_key in ("v1", "v2"):
            entry = result[engine_key]
            metrics = entry["metrics"]
            domestication = entry["domestication"]
            validation = entry["validation"]
            row = (
                engine_key,
                str(metrics["length_bp"]),
                f"{metrics['cai']:.3f}",
                f"{metrics['gc']:.1f}",
                str(metrics["polya_count"]),
                str(metrics["restriction_sites"]),
                str(metrics["splice_sites"]),
                str(domestication["removed_sites"]),
                "yes" if domestication["success"] else "no",
                "yes" if validation["valid"] else "no",
            )
            lines.append(" | ".join(row))
        lines.append("")

    return "\n".join(lines).rstrip()


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare v1 vs v2 outputs for protein FASTA")
    parser.add_argument("fasta", nargs="?", type=Path, help="Protein FASTA input")
    parser.add_argument("--input", dest="input_path", type=Path, help="Protein FASTA input")
    parser.add_argument("--out", dest="json_out", type=Path, help="Optional JSON output file")
    parser.add_argument("--profile", default="balanced", help="v2 profile (balanced/high_cai/gc_target)")
    parser.add_argument("--assembly", default="golden_gate", help="Restriction scan standard")
    args = parser.parse_args()

    fasta_path = args.input_path or args.fasta
    if fasta_path is None:
        raise SystemExit("Missing FASTA input. Provide positional or --input.")

    sequences = load_fasta(fasta_path)
    if not sequences:
        raise SystemExit("No FASTA records found.")

    try:
        profile = OptimizationProfile(args.profile.lower())
    except ValueError:
        profile = OptimizationProfile.BALANCED

    config = ComparisonConfig(profile=profile, assembly_standard=args.assembly)
    results = compare_sequences(sequences, config)

    print(render_table(results))

    payload = json.dumps(results, indent=2)
    if args.json_out:
        args.json_out.write_text(payload, encoding="utf-8")
    else:
        print("\nJSON:")
        print(payload)


if __name__ == "__main__":
    main()
