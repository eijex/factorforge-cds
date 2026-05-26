#requires -Version 5.1
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
Set-Location $RepoRoot
Write-Host "==> CWD: $(Get-Location)"

$FailedStep = $null

function Invoke-Step {
  param(
    [Parameter(Mandatory = $true)][string]$Name,
    [Parameter(Mandatory = $true)][scriptblock]$Action
  )
  $script:FailedStep = $Name
  Write-Host "==> $Name"
  & $Action
  if ($LASTEXITCODE -ne 0) {
    throw "$Name failed with exit code $LASTEXITCODE"
  }
}

try {
  $venvDir = ".venv-v2-panel"
  if (Test-Path $venvDir) {
    Remove-Item -Recurse -Force $venvDir
  }

  Invoke-Step "Create venv (.venv-v2-panel)" { python -m venv $venvDir }

  $activate = Join-Path $venvDir "Scripts\Activate.ps1"
  if (-not (Test-Path $activate)) {
    throw "Activate script not found: $activate"
  }
  Invoke-Step "Activate venv" { . $activate }

  $py = Join-Path $venvDir "Scripts\python.exe"

  Invoke-Step "Upgrade pip" { & $py -m pip install -U pip }

  $installFailed = $false
  try {
    Invoke-Step "Install (v2,dev) editable" { & $py -m pip install -e ".[v2,dev]" }
  } catch {
    $installFailed = $true
    Write-Host "Install (v2,dev) failed. Falling back to (dev)."
    Write-Host $_.Exception.Message
  }

  if ($installFailed) {
    Invoke-Step "Install (dev) editable" { & $py -m pip install -e ".[dev]" }
  }

  $panelScript = Join-Path $RepoRoot "scripts\audit\fetch_uniprot_panel.ps1"
  Invoke-Step "Download UniProt panel" { & $panelScript }

  $panelPath = Join-Path $RepoRoot "scripts\audit\inputs\uniprot_panel\panel_proteins.fasta"
  $profiles = @("balanced", "high_cai", "gc_target", "assembly_friendly", "ramp")

  New-Item -ItemType Directory -Force -Path "artifacts" | Out-Null

  $codonforgeExe = Join-Path $venvDir "Scripts\codonforge.exe"
  foreach ($profile in $profiles) {
    $outPath = "artifacts/v2.panel.$profile.fasta"
    if (Test-Path $codonforgeExe) {
      Invoke-Step "Optimize v2 ($profile)" { & $codonforgeExe optimize $panelPath --engine v2 --profile $profile --output $outPath }
    } else {
      Invoke-Step "Optimize v2 ($profile)" { & $py -m codonforge.cli.main optimize $panelPath --engine v2 --profile $profile --output $outPath }
    }
  }

  Invoke-Step "Audit v2 panel" { & $py scripts/audit/audit_v2_panel.py --input $panelPath --out artifacts/v2.panel.audit.json }

  if ($env:SKIP_PYTEST -eq "1") {
    Write-Host "Skipping pytest (SKIP_PYTEST=1)"
  } else {
    Invoke-Step "Pytest profile engine" { & $py -m pytest -q tests/engines/profile }
  }

  Write-Host "Outputs:"
  Write-Host "  $panelPath"
  Write-Host ("  " + (Join-Path $RepoRoot 'artifacts\v2.panel.balanced.fasta'))
  Write-Host ("  " + (Join-Path $RepoRoot 'artifacts\v2.panel.high_cai.fasta'))
  Write-Host ("  " + (Join-Path $RepoRoot 'artifacts\v2.panel.gc_target.fasta'))
  Write-Host ("  " + (Join-Path $RepoRoot 'artifacts\v2.panel.assembly_friendly.fasta'))
  Write-Host ("  " + (Join-Path $RepoRoot 'artifacts\v2.panel.ramp.fasta'))
  Write-Host ("  " + (Join-Path $RepoRoot 'artifacts\v2.panel.audit.json'))

  Write-Host "V2 PANEL ONE-SHOT: PASS"
  exit 0
} catch {
  Write-Host "V2 PANEL ONE-SHOT: FAIL"
  if ($FailedStep) {
    Write-Host "FAILED STEP: $FailedStep"
  }
  Write-Host $_.Exception.Message
  exit 1
}
