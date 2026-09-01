"""Run mSigAct signature assignment for all gene spectra.

mSigAct is an R package, so this script prepares the input data and calls the R script.

Usage:
  python 2_run_msigact.py
"""

from pathlib import Path
import subprocess
import sys
import json
import pandas as pd

HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "input"
RESULTS_DIR = HERE / "results"
WORK_DIR = HERE / "work" / "msigact"

BRANCHES = ("low_Ts", "high_Ts", "high_minus_low_Ts")
SPECTRUM_TYPES = ("Ts_only", "Ts_and_Tv")


def prepare_msigact_input() -> None:
    """Combine prepared branch matrices into a single input file for mSigAct."""
    branch_tables = []
    mutation_types = None

    for branch in BRANCHES:
        path = INPUT_DIR / f"{branch}_samples.txt"
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}; run 1_run_sigprofiler.py first")

        table = pd.read_csv(path, sep="\t", index_col=0)
        if mutation_types is None:
            mutation_types = table.index.tolist()
        elif mutation_types != table.index.tolist():
            raise ValueError("Branch matrices use different SBS96 row orders")

        # Rename columns to include branch prefix
        table.columns = [f"{branch}__{col}" for col in table.columns]
        branch_tables.append(table)

    # Combine all branches
    combined = pd.concat(branch_tables, axis=1)
    combined.index.name = "MutationType"

    # Save combined matrix
    output_path = INPUT_DIR / "samples_mSigAct.txt"
    combined.to_csv(output_path, sep="\t")
    print(f"✓ Prepared mSigAct input: {output_path}")
    print(f"  Dimensions: {combined.shape[0]} contexts × {combined.shape[1]} samples")


def run_r_script() -> bool:
    """Try to run the R script for mSigAct analysis."""
    r_script = HERE / "run_msigact.R"

    if not r_script.is_file():
        print(f"✗ R script not found: {r_script}")
        return False

    try:
        result = subprocess.run(
            ["Rscript", str(r_script)],
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            print(f"✗ R script failed with code {result.returncode}")
            print(result.stdout)
            print(result.stderr)
            return False

        print(result.stdout)
        return True

    except FileNotFoundError:
        print("✗ Rscript not found. Install R and add it to PATH.")
        return False
    except subprocess.TimeoutExpired:
        print("✗ R script timeout (>1 hour)")
        return False
    except Exception as e:
        print(f"✗ Error running R script: {e}")
        return False


def verify_results() -> bool:
    """Verify that mSigAct results were created."""
    results = [
        RESULTS_DIR / "msigact_activities.csv",
        RESULTS_DIR / "msigact_quality.csv",
        RESULTS_DIR / "msigact_empirical_priors.csv",
    ]

    all_exist = all(r.is_file() for r in results)
    if all_exist:
        print("\n✓ mSigAct results found:")
        for result in results:
            size = result.stat().st_size
            print(f"  - {result.name} ({size:,} bytes)")
        return True

    return False


def main() -> None:
    print("mSigAct Signature Assignment")
    print("=" * 60)

    print("\nStep 1: Preparing input...")
    try:
        prepare_msigact_input()
    except Exception as e:
        print(f"✗ Failed to prepare input: {e}")
        sys.exit(1)

    print("\nStep 2: Checking for existing results...")
    if verify_results():
        print("\n✓ Results already exist. Skipping R script execution.")
        print("\nTo re-run the analysis, delete the results files:")
        print(f"  rm {RESULTS_DIR / 'msigact*.csv'}")
        return

    print("\nStep 3: Running mSigAct (R)...")
    if not run_r_script():
        print("\n✗ mSigAct execution failed.")
        print("\nTo run mSigAct manually:")
        print(f"  Rscript {HERE / 'run_msigact.R'}")
        sys.exit(1)

    print("\nStep 4: Verifying results...")
    if verify_results():
        print("\n✓ mSigAct analysis complete!")
    else:
        print("\n✗ Results were not created properly.")
        sys.exit(1)


if __name__ == "__main__":
    main()
