#requires -Version 5.1
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$InputDir = Join-Path $RepoRoot "scripts\audit\inputs\uniprot_panel"
$AccessionsPath = Join-Path $InputDir "accessions.txt"
$MergedPath = Join-Path $InputDir "panel_proteins.fasta"

if (-not (Test-Path $AccessionsPath)) {
  throw "Accessions file not found: $AccessionsPath"
}

New-Item -ItemType Directory -Force -Path $InputDir | Out-Null

$accessions = @(
  Get-Content -Path $AccessionsPath |
    ForEach-Object { $_.Trim() } |
    Where-Object { $_ -and -not $_.StartsWith("#") }
)

if ($accessions.Count -eq 0) {
  throw "No accessions found in $AccessionsPath"
}

$downloaded = 0
foreach ($acc in $accessions) {
  $url = "https://rest.uniprot.org/uniprotkb/$acc.fasta"
  $outFile = Join-Path $InputDir "$acc.fasta"
  Invoke-WebRequest -Uri $url -OutFile $outFile
  $downloaded += 1
}

$builder = New-Object System.Text.StringBuilder
foreach ($acc in $accessions) {
  $src = Join-Path $InputDir "$acc.fasta"
  if (-not (Test-Path $src)) {
    throw "Missing downloaded FASTA: $src"
  }
  $content = Get-Content -Raw -Path $src
  [void]$builder.Append($content)
  if (-not $content.EndsWith("`n")) {
    [void]$builder.AppendLine()
  }
}

[System.IO.File]::WriteAllText($MergedPath, $builder.ToString(), [System.Text.UTF8Encoding]::new($false))

Write-Host "Downloaded FASTA files: $downloaded"
Write-Host "Merged panel: $MergedPath"
