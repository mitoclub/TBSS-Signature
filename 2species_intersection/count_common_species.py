"""Count species shared by every gene combination in the spectrum dataset.

The script uses the ready 12-component mutation spectra rather than raw sequence
metadata.  This matters because a species can have a sequence for a gene while
still lacking a mutation spectrum that can be used in the downstream comparison.

By default, inclusive intersections are calculated for ``Mammalia`` only.
"Inclusive" means that a species shared by four genes is also counted in each
relevant pair and triple.  This is the sample size needed when a particular
combination is analysed.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    PROJECT_ROOT / "1init_data" / "data" / "MutSpecVertebrates12.csv.gz"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "2species_intersection" / "data"
)
REQUIRED_COLUMNS = {"Gene", "Class", "Species", "Mut", "MutSpec"}
EXPECTED_MUTATIONS = frozenset(
    {
        "A>C",
        "A>G",
        "A>T",
        "C>A",
        "C>G",
        "C>T",
        "G>A",
        "G>C",
        "G>T",
        "T>A",
        "T>C",
        "T>G",
    }
)
DEFAULT_TAXONOMIC_CLASS = "Mammalia"


def load_gene_species_profiles(
    input_path: Path,
    taxonomic_class: str | None = DEFAULT_TAXONOMIC_CLASS,
) -> pd.DataFrame:
    """Return one row per gene/species profile in the requested scope.

    Gene and species identifiers are matched exactly as stored in the spectrum
    file.  We intentionally do not replace underscores, change case, or merge
    names because that could conflate distinct identifiers without a taxonomy
    reconciliation table.
    """

    spectra = pd.read_csv(
        input_path,
        usecols=lambda column: column in REQUIRED_COLUMNS,
    )

    missing_columns = REQUIRED_COLUMNS.difference(spectra.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing required columns in {input_path}: {missing}")

    if taxonomic_class is None:
        selected_spectra = spectra.copy()
    else:
        requested_class = taxonomic_class.strip()
        if not requested_class:
            raise ValueError("taxonomic_class must not be empty")
        class_mask = (
            spectra["Class"].astype("string").str.strip().str.casefold()
            == requested_class.casefold()
        )
        selected_spectra = spectra.loc[class_mask].copy()
        if selected_spectra.empty:
            available_classes = sorted(
                spectra["Class"].dropna().astype(str).unique()
            )
            raise ValueError(
                f"No rows found for class {taxonomic_class!r}. "
                f"Available classes: {available_classes}"
            )

    if selected_spectra.empty:
        raise ValueError(f"No spectrum rows found in {input_path}")
    if selected_spectra[["Gene", "Species", "Mut"]].isna().any().any():
        raise ValueError("Gene, Species, and Mut identifiers must not be missing")
    if selected_spectra.duplicated(["Gene", "Species", "Mut"]).any():
        raise ValueError("Duplicate Gene/Species/Mut rows were found")

    mutation_sets = selected_spectra.groupby(["Gene", "Species"])["Mut"].agg(
        frozenset
    )
    invalid_mutation_sets = mutation_sets[mutation_sets.ne(EXPECTED_MUTATIONS)]
    if not invalid_mutation_sets.empty:
        raise ValueError(
            f"Found {len(invalid_mutation_sets)} Gene/Species profiles without "
            "the exact 12 substitution categories"
        )

    mutspec = pd.to_numeric(selected_spectra["MutSpec"], errors="coerce")
    if mutspec.isna().any() or not np.isfinite(mutspec).all():
        raise ValueError("MutSpec values must be numeric, complete, and finite")
    if (mutspec < 0).any():
        raise ValueError("MutSpec values must be non-negative")

    profile_sums = mutspec.groupby(
        [selected_spectra["Gene"], selected_spectra["Species"]]
    ).sum()
    if not np.allclose(profile_sums.to_numpy(), 1.0, rtol=0, atol=1e-8):
        raise ValueError("MutSpec must sum to one in every Gene/Species profile")

    # A spectrum contains one row per mutation type.  Dropping duplicates here
    # reduces it to the presence/absence table needed for set intersections.
    # Class is retained so the Mammalia-only scope remains auditable.
    return (
        selected_spectra[["Gene", "Species", "Class"]]
        .drop_duplicates()
        .sort_values(["Gene", "Species"])
    )


def build_species_sets(profiles: pd.DataFrame) -> dict[str, set[str]]:
    """Map every gene to the species for which a spectrum is available."""

    return {
        gene: set(group["Species"])
        for gene, group in profiles.groupby("Gene", sort=True)
    }


def count_intersections(
    species_by_gene: dict[str, set[str]],
    max_combination_size: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build inclusive intersection counts and a long table of their members."""

    genes = sorted(species_by_gene)
    largest_size = min(max_combination_size, len(genes))
    count_records: list[dict[str, object]] = []
    member_records: list[dict[str, object]] = []

    for size in range(2, largest_size + 1):
        for gene_combination in combinations(genes, size):
            shared_species = set.intersection(
                *(species_by_gene[gene] for gene in gene_combination)
            )
            union_species = set.union(
                *(species_by_gene[gene] for gene in gene_combination)
            )
            combination_label = " + ".join(gene_combination)

            count_records.append(
                {
                    "combination_size": size,
                    "genes": combination_label,
                    "n_common_species": len(shared_species),
                    "n_union_species": len(union_species),
                    "jaccard_index": (
                        len(shared_species) / len(union_species)
                        if union_species
                        else float("nan")
                    ),
                }
            )

            # Keeping the member table makes every reported count auditable and
            # lets downstream analyses select exactly the same matched species.
            member_records.extend(
                {
                    "combination_size": size,
                    "genes": combination_label,
                    "species": species,
                }
                for species in sorted(shared_species)
            )

    counts = pd.DataFrame.from_records(
        count_records,
        columns=[
            "combination_size",
            "genes",
            "n_common_species",
            "n_union_species",
            "jaccard_index",
        ],
    )
    members = pd.DataFrame.from_records(
        member_records,
        columns=["combination_size", "genes", "species"],
    )
    return counts, members


def make_gene_summary(species_by_gene: dict[str, set[str]]) -> pd.DataFrame:
    """Summarise how many usable species are available for each gene."""

    return pd.DataFrame(
        {
            "gene": sorted(species_by_gene),
            "n_species": [
                len(species_by_gene[gene]) for gene in sorted(species_by_gene)
            ],
        }
    ).sort_values(["n_species", "gene"], ascending=[False, True])


def print_results(
    gene_summary: pd.DataFrame,
    intersection_counts: pd.DataFrame,
    taxonomic_scope: str,
) -> None:
    """Print compact human-readable tables to the terminal."""

    print(f"\nSpecies with usable spectra in scope: {taxonomic_scope}")
    print(gene_summary.to_string(index=False))

    if intersection_counts.empty:
        print("\nNo multi-gene intersections are available.")
        return

    for size, table in intersection_counts.groupby("combination_size", sort=True):
        print(f"\nInclusive {size}-gene intersections")
        print(
            table.loc[:, ["genes", "n_common_species"]]
            .sort_values(["n_common_species", "genes"], ascending=[False, True])
            .to_string(index=False)
        )


def parse_args(args: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"12-component spectrum CSV (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--taxonomic-class",
        default=DEFAULT_TAXONOMIC_CLASS,
        help=(
            "Exact value from the Class column "
            f"(default: {DEFAULT_TAXONOMIC_CLASS})."
        ),
    )
    parser.add_argument(
        "--max-combination-size",
        type=int,
        default=4,
        help="Largest gene combination to report (default: 4)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for result tables (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args(args)


def main(args: Iterable[str] | None = None) -> None:
    options = parse_args(args)
    if options.max_combination_size < 2:
        raise ValueError("--max-combination-size must be at least 2")

    input_path = options.input.expanduser().resolve()
    output_dir = options.output_dir.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Spectrum file not found: {input_path}")

    profiles = load_gene_species_profiles(input_path, options.taxonomic_class)
    species_by_gene = build_species_sets(profiles)
    if len(species_by_gene) < 2:
        raise ValueError(
            "At least two genes are required to calculate an intersection"
        )
    gene_summary = make_gene_summary(species_by_gene)
    intersection_counts, intersection_members = count_intersections(
        species_by_gene,
        max_combination_size=options.max_combination_size,
    )
    species_classes = (
        profiles[["Species", "Class"]]
        .drop_duplicates()
        .rename(columns={"Species": "species", "Class": "taxonomic_class"})
    )
    if species_classes["species"].duplicated().any():
        raise ValueError("A species identifier occurs in more than one class")
    intersection_members = intersection_members.merge(
        species_classes,
        on="species",
        how="left",
        validate="many_to_one",
    )

    try:
        # Forward slashes so the saved tables are identical on every platform.
        source_label = input_path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        source_label = input_path.as_posix()
    scope_label = options.taxonomic_class.strip()
    for table in (gene_summary, intersection_counts, intersection_members):
        table.insert(0, "source_file", source_label)
        table.insert(0, "taxonomic_scope", scope_label)

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_summary.to_csv(output_dir / "gene_species_counts.csv", index=False)
    intersection_counts.to_csv(
        output_dir / "species_intersection_counts.csv", index=False
    )
    intersection_members.to_csv(
        output_dir / "species_intersection_members.csv", index=False
    )

    print_results(gene_summary, intersection_counts, scope_label)
    print(f"\nTables written to: {output_dir}")


if __name__ == "__main__":
    main()
