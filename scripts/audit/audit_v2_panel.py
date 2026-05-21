#!/usr/bin/env python3
"""Audit FactorForge v2 panel outputs and emit JSON summary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from factorforge.engines.registry import EngineRegistry


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    seq_id: str | None = None
    seq_lines: list[str] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id is not None:
                records.append((seq_id, "".join(seq_lines)))
            header = line[1:].strip()
            seq_id = header.split()[0] if header else f"seq{len(records) + 1}"
            seq_lines = []
        else:
            seq_lines.append(line)

    if seq_id is not None:
        records.append((seq_id, "".join(seq_lines)))

    return records


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for base in seq if base.upper() in {"G", "C"})
    return 100.0 * gc / len(seq)


def metrics_to_dict(metrics: Any) -> dict[str, Any]:
    if metrics is None:
        return {}
    if isinstance(metrics, dict):
        return dict(metrics)
    if hasattr(metrics, "dict") and callable(metrics.dict):
        try:
            return dict(metrics.dict())
        except Exception:
            return {"raw": str(metrics)}
    try:
        return dict(metrics)
    except Exception:
        return {"raw": str(metrics)}


def summarize(stats: dict[str, dict[str, list[float]]]) -> None:
    print("Audit summary")
    for profile, values in stats.items():
        lengths = values["lengths"]
        gcs = values["gcs"]
        if not lengths:
            print(f"  {profile}: count=0")
            continue
        count = len(lengths)
        mean_gc = sum(gcs) / count
        min_len = int(min(lengths))
        max_len = int(max(lengths))
        print(
            f"  {profile}: count={count} mean_gc={mean_gc:.2f} "
            f"min_len_bp={min_len} max_len_bp={max_len}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit v2 panel outputs")
    parser.add_argument("--input", required=True, help="Input protein FASTA")
    parser.add_argument("--out", required=True, help="Output JSON path")
    parser.add_argument(
        "--profiles",
        default="balanced,high_cai,gc_target,assembly_friendly,ramp",
        help="Comma-separated profiles",
    )
    parser.add_argument("--engine", default="v2", help="Engine name")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.out)

    if not input_path.exists():
        raise SystemExit(f"Input FASTA not found: {input_path}")

    profiles = [p.strip() for p in args.profiles.split(",") if p.strip()]
    if not profiles:
        raise SystemExit("No profiles specified")

    records = parse_fasta(input_path)
    if not records:
        raise SystemExit("No FASTA records found")

    engine = EngineRegistry.get(args.engine)

    stats: dict[str, dict[str, list[float]]] = {
        profile: {"lengths": [], "gcs": []} for profile in profiles
    }

    results: list[dict[str, Any]] = []
    for seq_id, protein_seq in records:
        entry: dict[str, Any] = {
            "id": seq_id,
            "protein_len": len(protein_seq),
            "outputs": {},
        }
        for profile in profiles:
            try:
                result = engine.optimize(protein_seq, profile=profile)
                dna = result.sequence
                metrics = metrics_to_dict(result.metrics)
                metrics["length_bp"] = len(dna) if dna else 0
                metrics["gc_percent"] = gc_percent(dna)
                entry["outputs"][profile] = {"dna": dna, "metrics": metrics}
                stats[profile]["lengths"].append(float(metrics["length_bp"]))
                stats[profile]["gcs"].append(float(metrics["gc_percent"]))
            except Exception as exc:
                entry["outputs"][profile] = {"error": str(exc)}
        results.append(entry)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    summarize(stats)
    print(f"Wrote audit JSON: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
