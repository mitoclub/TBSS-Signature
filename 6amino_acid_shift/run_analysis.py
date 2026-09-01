"""Extract and QC the mammalian MIDORI2 data used by stage 6."""

from __future__ import annotations

import argparse
from pathlib import Path

from midori_analysis import run_midori_analysis


def repository_root() -> Path:
    here = Path.cwd().resolve()
    for candidate in (here, *here.parents):
        if (candidate / "6amino_acid_shift" / "data" / "midori").is_dir():
            return candidate
    raise FileNotFoundError("Run from the repository root or 6amino_acid_shift directory")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--class-filter",
        default="Mammalia",
        choices=("Mammalia",),
        help="MIDORI taxonomic class (currently Mammalia only)",
    )
    args = parser.parse_args()

    root = repository_root()
    result = run_midori_analysis(
        root / "6amino_acid_shift",
        class_filter=args.class_filter,
    )
    summary = result["summary"]
    print()
    print(f"Taxonomic scope: {result['taxonomic_scope']}")
    print(
        "MIDORI2 taxa represented by at least one target gene: "
        f"{summary['n_taxa_any_target_gene']}"
    )
    print(
        "Taxa with nucleotide records for all target genes: "
        f"{summary['n_taxa_all_target_nucleotide_genes']}"
    )
    print(
        "Taxa in the final paired QC cohort: "
        f"{summary['n_taxa_final_paired_qc']}"
    )
    print(f"Retrieval backend: {result['retrieval_backend']}")


if __name__ == "__main__":
    main()
