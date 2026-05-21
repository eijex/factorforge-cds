#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INPUT_DIR="$ROOT_DIR/scripts/audit/inputs/uniprot_panel"
ACCESSIONS_FILE="$INPUT_DIR/accessions.txt"
MERGED_FILE="$INPUT_DIR/panel_proteins.fasta"

if [ ! -f "$ACCESSIONS_FILE" ]; then
  echo "Accessions file not found: $ACCESSIONS_FILE" >&2
  exit 1
fi

mkdir -p "$INPUT_DIR"

if command -v curl >/dev/null 2>&1; then
  fetch() {
    local url="$1"
    local out="$2"
    curl -fsSL "$url" -o "$out"
  }
elif command -v wget >/dev/null 2>&1; then
  fetch() {
    local url="$1"
    local out="$2"
    wget -q -O "$out" "$url"
  }
else
  echo "Neither curl nor wget is available." >&2
  exit 1
fi

ACCESSIONS=()
while IFS= read -r line; do
  line="$(printf '%s' "$line" | sed -e 's/\r$//' -e 's/#.*$//' -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
  if [ -n "$line" ]; then
    ACCESSIONS+=("$line")
  fi
done < "$ACCESSIONS_FILE"

if [ ${#ACCESSIONS[@]} -eq 0 ]; then
  echo "No accessions found in $ACCESSIONS_FILE" >&2
  exit 1
fi

count=0
for acc in "${ACCESSIONS[@]}"; do
  url="https://rest.uniprot.org/uniprotkb/${acc}.fasta"
  out="$INPUT_DIR/${acc}.fasta"
  fetch "$url" "$out"
  count=$((count + 1))
done

: > "$MERGED_FILE"
for acc in "${ACCESSIONS[@]}"; do
  src="$INPUT_DIR/${acc}.fasta"
  if [ ! -f "$src" ]; then
    echo "Missing downloaded FASTA: $src" >&2
    exit 1
  fi
  cat "$src" >> "$MERGED_FILE"
  printf '\n' >> "$MERGED_FILE"
done

echo "Downloaded FASTA files: $count"
echo "Merged panel: $MERGED_FILE"
