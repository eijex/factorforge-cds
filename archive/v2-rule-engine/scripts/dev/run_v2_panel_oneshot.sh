#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

VENV_DIR=".venv-v2-panel"
PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  PYTHON_BIN="python"
fi

FAIL_CMD=""
FAIL_OUTPUT=""
FAIL_SUGGESTION=""
LAST_OUTPUT=""

run_step() {
  local description="$1"
  local suggestion="$2"
  shift 2
  local cmd=("$@")

  printf "==> %s\n" "$description"
  set +e
  local output
  output="$(${cmd[@]} 2>&1)"
  local status=$?
  set -e

  LAST_OUTPUT="$output"

  if [ $status -ne 0 ]; then
    FAIL_CMD="${cmd[*]}"
    FAIL_OUTPUT="$output"
    FAIL_SUGGESTION="$suggestion"
    return $status
  fi

  return 0
}

fail_and_exit() {
  echo "V2 PANEL ONE-SHOT: FAIL"
  echo "FAILED COMMAND: $FAIL_CMD"
  echo "OUTPUT:"
  echo "$FAIL_OUTPUT"
  echo "SUGGESTION: $FAIL_SUGGESTION"
  exit 1
}

if [ -d "$VENV_DIR" ]; then
  rm -rf "$VENV_DIR"
fi

if ! run_step "Create fresh v2 panel venv" "Ensure Python 3 with venv support is installed." "$PYTHON_BIN" -m venv "$VENV_DIR"; then
  fail_and_exit
fi

# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"

if ! run_step "Upgrade pip" "Check network access to PyPI." python -m pip install -U pip; then
  fail_and_exit
fi

if ! run_step "Install (v2,dev) editable" "Check network access and extras in pyproject.toml." python -m pip install -e ".[v2,dev]"; then
  echo "Install (v2,dev) failed. Falling back to (dev)."
  echo "$LAST_OUTPUT"
  if ! run_step "Install (dev) editable" "Check network access and extras in pyproject.toml." python -m pip install -e ".[dev]"; then
    fail_and_exit
  fi
fi

export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8

if ! run_step "Download UniProt panel" "Check network access to rest.uniprot.org." bash scripts/audit/fetch_uniprot_panel.sh; then
  fail_and_exit
fi

panel_path="scripts/audit/inputs/uniprot_panel/panel_proteins.fasta"
mkdir -p artifacts

if command -v codonforge >/dev/null 2>&1; then
  run_optimize() {
    codonforge optimize "$1" --engine v2 --profile "$2" --output "$3"
  }
else
  run_optimize() {
    python -m codonforge.cli.main optimize "$1" --engine v2 --profile "$2" --output "$3"
  }
fi

for profile in balanced high_cai gc_target assembly_friendly ramp; do
  out_path="artifacts/v2.panel.${profile}.fasta"
  if ! run_step "Optimize v2 (${profile})" "Verify profile engine installs and the panel FASTA is present." run_optimize "$panel_path" "$profile" "$out_path"; then
    fail_and_exit
  fi
  if [ -n "$LAST_OUTPUT" ]; then
    printf "%s\n" "$LAST_OUTPUT"
  fi
done

if ! run_step "Audit v2 panel" "Verify the audit script and panel FASTA are present." python scripts/audit/audit_v2_panel.py --input "$panel_path" --out artifacts/v2.panel.audit.json; then
  fail_and_exit
fi
if [ -n "$LAST_OUTPUT" ]; then
  printf "%s\n" "$LAST_OUTPUT"
fi

if [ "${SKIP_PYTEST:-}" = "1" ]; then
  echo "Skipping pytest (SKIP_PYTEST=1)"
else
  if ! run_step "Pytest profile engine" "Inspect failing tests or ensure dependencies are installed." python -m pytest -q tests/engines/profile; then
    fail_and_exit
  fi
fi

echo "Outputs:"
echo "  $ROOT_DIR/$panel_path"
echo "  $ROOT_DIR/artifacts/v2.panel.balanced.fasta"
echo "  $ROOT_DIR/artifacts/v2.panel.high_cai.fasta"
echo "  $ROOT_DIR/artifacts/v2.panel.gc_target.fasta"
echo "  $ROOT_DIR/artifacts/v2.panel.assembly_friendly.fasta"
echo "  $ROOT_DIR/artifacts/v2.panel.ramp.fasta"
echo "  $ROOT_DIR/artifacts/v2.panel.audit.json"

echo "V2 PANEL ONE-SHOT: PASS"
