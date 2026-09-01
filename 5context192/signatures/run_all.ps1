$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = Join-Path $scriptDir ".venv\Scripts\python.exe"
if (-not (Test-Path -LiteralPath $python)) {
    throw "Missing local Python environment. Create .venv and install requirements.txt."
}

$rCandidates = @(
    "C:\Program Files\R\R-4.6.1\bin\Rscript.exe",
    "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
)
$rscript = $rCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
if ($null -eq $rscript) {
    $rCommand = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($null -eq $rCommand) {
        throw "Rscript was not found."
    }
    $rscript = $rCommand.Source
}

& $python (Join-Path $scriptDir "prepare_inputs.py")
& $python (Join-Path $scriptDir "run_sigprofiler_assignment.py")
& $rscript (Join-Path $scriptDir "run_msigact.R")
& $python (Join-Path $scriptDir "summarize_results.py")
& $python (Join-Path $scriptDir "verify_signature_analysis.py")
