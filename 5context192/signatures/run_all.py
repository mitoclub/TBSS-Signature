"""Run the whole COSMIC signature substage on any platform.

    python run_all.py

Steps, in order:

1. ``1_run_sigprofiler.py`` prepares the SBS96 inputs and runs
   SigProfilerAssignment;
2. ``2_run_msigact.py`` builds the mSigAct input and calls ``run_msigact.R``
   through ``Rscript``.

Both steps run with the interpreter that started this script, so activating the
substage virtual environment first is enough.  ``Rscript`` must be on PATH;
install the R packages once with ``Rscript install_r_dependencies.R``.
"""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STEPS = ("1_run_sigprofiler.py", "2_run_msigact.py")


def main() -> int:
    if shutil.which("Rscript") is None:
        print(
            "Rscript was not found on PATH. Install R, then run "
            "`Rscript install_r_dependencies.R` before this script.",
            file=sys.stderr,
        )
        return 1

    for step in STEPS:
        script = HERE / step
        if not script.is_file():
            print(f"Missing step script: {script}", file=sys.stderr)
            return 1
        print(f"\n=== {step} ===", flush=True)
        result = subprocess.run([sys.executable, str(script)], cwd=HERE)
        if result.returncode != 0:
            print(f"{step} failed with exit code {result.returncode}", file=sys.stderr)
            return result.returncode

    print("\nSignature substage finished. Results are in results/.")
    print("Open SignatureAnalysis.ipynb for the report and figures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
