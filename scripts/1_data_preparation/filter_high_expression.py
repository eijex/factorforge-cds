"""Filter high-expression N. benthamiana proteins from Salmon TPM output."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def filter_by_tpm(
    quant_path: str,
    output_path: str,
    percentile: float = 75.0,
) -> None:
    """Write IDs with TPM values at or above the selected percentile."""
    rows: list[tuple[str, float]] = []
    with open(quant_path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tpm = float(row["TPM"])
            if tpm > 0:
                rows.append((row["Name"], tpm))

    if not rows:
        raise ValueError(f"No TPM values found in {quant_path}")

    tpms = sorted(tpm for _, tpm in rows)
    cutoff_index = min(int(len(tpms) * percentile / 100), len(tpms) - 1)
    cutoff = tpms[cutoff_index]
    high_expression = [name for name, tpm in rows if tpm >= cutoff]

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(high_expression), encoding="utf-8")

    print(f"TPM cutoff ({percentile}th pct): {cutoff:.2f}")
    print(f"High-expression proteins: {len(high_expression)}")
    print(f"Saved to: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quant", default="data/raw/grosse_holz_quant.sf")
    parser.add_argument("--output", default="data/raw/high_expression_proteins.txt")
    parser.add_argument("--percentile", type=float, default=75.0)
    args = parser.parse_args()
    filter_by_tpm(args.quant, args.output, args.percentile)


if __name__ == "__main__":
    main()
