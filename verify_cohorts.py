"""Verify the matched-cohort and gene-order contract of stages 2-5.

Stages 2-5 must compare genes only on species that have a complete, valid
profile for *every* selected gene, and every stage must place genes in the one
canonical mitochondrial order defined in ``mtdna.py``.  This script re-derives
both properties from the source tables and compares them with the tables saved
in each stage, so a silent change in cohort composition cannot pass unnoticed.

Run it from the repository root or any stage directory::

    python verify_cohorts.py

Stage 6 has its own contract check in ``6amino_acid_shift/verify_midori.py``.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

_REPOSITORY_ROOT = Path(__file__).resolve().parent
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from mtdna import CANONICAL_ORDER, canonical_order, gene_metadata_frame

CLASS_FILTER = "Mammalia"
SPECTRUM_12 = Path("1init_data/data/MutSpecVertebrates12.csv.gz")
SPECTRUM_192 = Path("1init_data/data/MutSpecVertebrates192.csv.gz")

# stage directory -> (source table, genes, saved cohort file, components per profile)
COHORTS = {
    "3compare_t_genes": (
        SPECTRUM_12,
        ("CO1", "Cytb"),
        Path("3compare_t_genes/data/common_species_CO1_vs_Cytb.csv"),
        12,
    ),
    "4tsss_gradient": (
        SPECTRUM_12,
        ("CO1", "CO3", "Cytb"),
        Path("4tsss_gradient/data/matched_species_CO1_CO3_Cytb.csv"),
        12,
    ),
    "5context192": (
        SPECTRUM_192,
        ("CO1", "CO3", "Cytb"),
        Path("5context192/data/matched_species_CO1_CO3_Cytb.csv"),
        192,
    ),
}

PROXY_TABLES = (
    Path("4tsss_gradient/data/gene_tsss_proxy.csv"),
    Path("5context192/data/gene_tsss_proxy.csv"),
)


def _mammals(root: Path, source: Path) -> pd.DataFrame:
    table = pd.read_csv(root / source)
    return table.loc[table["Class"].eq(CLASS_FILTER)]


def _intersection(mammals: pd.DataFrame, genes: tuple[str, ...]) -> set[str]:
    sets = [
        set(mammals.loc[mammals["Gene"].eq(gene), "Species"])
        for gene in genes
    ]
    return set.intersection(*sets)


def check_gene_order() -> None:
    assert canonical_order(reversed(CANONICAL_ORDER)) == CANONICAL_ORDER
    for genes in (("Cytb", "CO1"), ("Cytb", "CO3", "CO1"), ("ND2", "CO1")):
        assert canonical_order(genes) == canonical_order(sorted(genes)), genes
    print("Gene order: canonical rCRS order is stable for any input order.")


def check_stage_two(root: Path) -> None:
    counts = pd.read_csv(root / "2species_intersection/data/species_intersection_counts.csv")
    mammals = _mammals(root, SPECTRUM_12)
    for row in counts.itertuples(index=False):
        genes = tuple(part.strip() for part in row.genes.split("+"))
        assert canonical_order(genes), genes
        recomputed = len(_intersection(mammals, genes))
        assert recomputed == row.n_common_species, (
            f"{row.genes}: saved {row.n_common_species}, recomputed {recomputed}"
        )
    print(f"Stage 2: all {len(counts)} saved intersections match the source table.")


def check_matched_cohorts(root: Path) -> None:
    for stage, (source, genes, cohort_file, components) in COHORTS.items():
        ordered = canonical_order(genes)
        assert ordered == genes, f"{stage} lists genes out of canonical order: {genes}"
        mammals = _mammals(root, source)
        expected = _intersection(mammals, ordered)
        saved = set(pd.read_csv(root / cohort_file)["Species"])
        assert saved == expected, (
            f"{stage}: saved cohort differs from the exact all-gene intersection "
            f"(saved {len(saved)}, expected {len(expected)})"
        )
        matched = mammals.loc[
            mammals["Species"].isin(saved) & mammals["Gene"].isin(ordered)
        ]
        profile_sizes = matched.groupby(["Species", "Gene"]).size()
        assert len(profile_sizes) == len(saved) * len(ordered), (
            f"{stage}: expected one profile per species and gene"
        )
        assert profile_sizes.eq(components).all(), (
            f"{stage}: every profile must contain {components} components"
        )
        print(
            f"Stage {stage[0]}: {len(saved):>3} species x {len(ordered)} genes "
            f"= exact intersection, {components} components per profile."
        )


def check_positional_metadata(root: Path) -> None:
    for table_path in PROXY_TABLES:
        saved = pd.read_csv(root / table_path)
        rebuilt = gene_metadata_frame(
            saved["Gene"], plot_label="plot_label" in saved.columns
        )
        assert list(saved["Gene"]) == list(rebuilt["Gene"]), table_path
        pd.testing.assert_frame_equal(
            saved.reset_index(drop=True),
            rebuilt[saved.columns].reset_index(drop=True),
            check_exact=False,
            rtol=1e-12,
        )
    print("Stages 4-5: saved DssH metadata matches mtdna.gene_metadata_frame.")


def main() -> None:
    here = Path.cwd().resolve()
    root = next(
        (path for path in (here, *here.parents) if (path / SPECTRUM_12).is_file()),
        None,
    )
    if root is None:
        raise FileNotFoundError("Run from the repository root or a stage directory")
    check_gene_order()
    check_stage_two(root)
    check_matched_cohorts(root)
    check_positional_metadata(root)
    print("Cohort verification passed.")


if __name__ == "__main__":
    main()
