"""Command-line entry point for the stage-6 MIDORI2 sequence/AA analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from midori_analysis import run_midori_analysis


def repository_root() -> Path:
    here = Path.cwd().resolve()
    for candidate in (here, *here.parents):
        if (candidate / "1init_data" / "data" / "MutSpecVertebrates12.csv.gz").exists():
            return candidate
    raise FileNotFoundError("Run from the repository root or 6amino_acid_shift directory")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    root = repository_root()
    common_arguments = {
        "mutation_source": root / "1init_data" / "data" / "MutSpecVertebrates12.csv.gz",
        "stage_dir": root / "6amino_acid_shift",
        "class_filter": "Mammalia",
        "n_boot": args.bootstrap,
    }
    result = run_midori_analysis(**common_arguments)
    attrition = result["attrition"].set_index("stage")["n_species"]
    n_mut = int(attrition["mutation-spectrum cohort"])
    n_seq = int(attrition["all 3 genes found"])
    n_final = int(attrition["final analysis cohort"])
    print()
    print(
        f"Mutation-spectrum data were available for all three genes in {n_mut} species.\n"
        f"MIDORI2 sequences for COX1, COX3 and CytB were successfully obtained for {n_seq} species.\n"
        f"After sequence QC, {n_final} species remained in the final matched cohort."
    )
    print(f"Retrieval backend: {result['retrieval_backend']}")


if __name__ == "__main__":
    main()
