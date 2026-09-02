"""Canonical mitochondrial gene order, coordinates, and names for every stage.

Mammalian mitochondrial gene order is conserved, so the project fixes one
canonical order here and derives it nowhere else.  Every stage imports this
module instead of repeating gene lists, display names, or rCRS coordinates:

* stages 2-5 analyse the ``MutSpec`` tables, which contain ``CO1``, ``CO3``,
  ``Cytb``, and ``ND2``;
* stage 6 analyses MIDORI2 records for the ten heavy-strand Major Arc genes.

Both naming systems use the same keys, so one table serves both.  Coordinates
are human rCRS (NC_012920.1) positions applied as a common positional scale;
they are not species-specific measurements.

Import it from any stage directory or the repository root::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(REPOSITORY_ROOT))
    from mtdna import MAJOR_ARC_HEAVY_STRAND, canonical_order, gene_metadata_frame
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

MTDNA_LENGTH = 16_569
"""Length of the human rCRS mitochondrial genome in base pairs."""

ORIL_START = 5_730
"""First rCRS position of the light-strand origin of replication."""


@dataclass(frozen=True)
class MitochondrialGene:
    """One protein-coding gene in the human rCRS reference."""

    key: str
    display: str
    start: int
    end: int
    strand: str

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2

    @property
    def dssh_proxy(self) -> float:
        """Positional proxy for the time spent single-stranded.

        ``2 * (midpoint - OriL) / genome length``.  The scale is linear and is
        not wrapped around the circular genome, so it is meaningful only for
        genes downstream of ``ORIL_START``.  Genes upstream of OriL (ND1, ND2)
        receive negative values, and CytB exceeds one; the number is an ordered
        positional coordinate, not a fraction of a replication cycle.
        """

        return 2 * (self.midpoint - ORIL_START) / MTDNA_LENGTH

    @property
    def is_downstream_of_oril(self) -> bool:
        return self.start >= ORIL_START


MTDNA_GENES: tuple[MitochondrialGene, ...] = (
    MitochondrialGene("ND1", "ND1", 3_307, 4_262, "H"),
    MitochondrialGene("ND2", "ND2", 4_470, 5_511, "H"),
    MitochondrialGene("CO1", "COX1", 5_904, 7_445, "H"),
    MitochondrialGene("CO2", "COX2", 7_586, 8_269, "H"),
    MitochondrialGene("A8", "ATP8", 8_366, 8_572, "H"),
    MitochondrialGene("A6", "ATP6", 8_527, 9_207, "H"),
    MitochondrialGene("CO3", "COX3", 9_207, 9_990, "H"),
    MitochondrialGene("ND3", "ND3", 10_059, 10_404, "H"),
    MitochondrialGene("ND4L", "ND4L", 10_470, 10_766, "H"),
    MitochondrialGene("ND4", "ND4", 10_760, 12_137, "H"),
    MitochondrialGene("ND5", "ND5", 12_337, 14_148, "H"),
    MitochondrialGene("ND6", "ND6", 14_149, 14_673, "L"),
    MitochondrialGene("Cytb", "CytB", 14_747, 15_887, "H"),
)
"""The 13 protein-coding genes in ascending rCRS start order."""

GENES_BY_KEY = {gene.key: gene for gene in MTDNA_GENES}

CANONICAL_ORDER: tuple[str, ...] = tuple(gene.key for gene in MTDNA_GENES)

DISPLAY_NAMES: dict[str, str] = {gene.key: gene.display for gene in MTDNA_GENES}

MAJOR_ARC_HEAVY_STRAND: tuple[str, ...] = tuple(
    gene.key
    for gene in MTDNA_GENES
    if gene.strand == "H" and gene.is_downstream_of_oril
)
"""COX1 to CytB: the ten heavy-strand genes analysed in stage 6.

ND6 lies in the same physical interval but is encoded on the light strand.
ND1 and ND2 lie upstream of OriL and are outside this Major Arc definition.
"""

# Ten-gene Major Arc palette used by stage 6.  Stages 3-5 keep their own
# three-gene palette in ``mutation_comparison.py`` so their published figures
# do not change colour.
MAJOR_ARC_GENE_COLORS: dict[str, str] = {
    "ND1": "#B07AA1",
    "ND2": "#CC79A7",
    "CO1": "#0072B2",
    "CO2": "#348ABD",
    "A8": "#E69F00",
    "A6": "#F0B44D",
    "CO3": "#56B4E9",
    "ND3": "#8E6CBE",
    "ND4L": "#9B7EBD",
    "ND4": "#7A5195",
    "ND5": "#5F4B8B",
    "ND6": "#D55E00",
    "Cytb": "#009E73",
}


def gene(key: str) -> MitochondrialGene:
    """Return one canonical gene record, with a helpful error for typos."""

    try:
        return GENES_BY_KEY[key]
    except KeyError:
        raise KeyError(
            f"Unknown mitochondrial gene {key!r}. "
            f"Known keys: {', '.join(CANONICAL_ORDER)}"
        ) from None


def canonical_order(genes: Iterable[str]) -> tuple[str, ...]:
    """Sort gene keys into canonical rCRS order.

    Any selection of genes is reordered the same way regardless of how it was
    typed, so figures, tables, and ordered trend tests cannot silently depend on
    the order used in a notebook.  Duplicates and unknown keys are rejected.
    """

    selected = tuple(genes)
    if len(set(selected)) != len(selected):
        raise ValueError(f"genes must be unique: {selected}")
    unknown = [key for key in selected if key not in GENES_BY_KEY]
    if unknown:
        raise KeyError(
            f"Unknown mitochondrial genes: {unknown}. "
            f"Known keys: {', '.join(CANONICAL_ORDER)}"
        )
    return tuple(sorted(selected, key=CANONICAL_ORDER.index))


def display_names(genes: Iterable[str]) -> tuple[str, ...]:
    """Return publication display names in the order supplied."""

    return tuple(gene(key).display for key in genes)


def dssh_proxy(key: str) -> float:
    """Positional single-stranded-duration proxy for one gene."""

    return gene(key).dssh_proxy


def arc_positions(genes: Sequence[str]) -> tuple[str, ...]:
    """Label the first, last, and intermediate genes of an ordered selection."""

    if len(genes) < 2:
        raise ValueError("arc positions need at least two genes")
    labels = ["middle"] * len(genes)
    labels[0] = "beginning"
    labels[-1] = "end"
    return tuple(labels)


def gene_metadata_frame(genes: Iterable[str], *, plot_label: bool = False):
    """Build the canonical positional metadata table used by stages 4 and 5.

    The genes are always returned in canonical rCRS order.  ``plot_label`` adds
    the two-line axis label used by the stage-4 figures.
    """

    import pandas as pd

    ordered = canonical_order(genes)
    records = [gene(key) for key in ordered]
    frame = pd.DataFrame(
        {
            "Gene": [record.key for record in records],
            "display_gene": [record.display for record in records],
            "arc_position": list(arc_positions(ordered)),
            "rCRS_start": [record.start for record in records],
            "rCRS_end": [record.end for record in records],
        }
    )
    frame["rCRS_midpoint"] = (frame["rCRS_start"] + frame["rCRS_end"]) / 2
    frame["dssh_proxy"] = 2 * (frame["rCRS_midpoint"] - ORIL_START) / MTDNA_LENGTH
    if plot_label:
        frame["plot_label"] = (
            frame["display_gene"]
            + "\nDssH="
            + frame["dssh_proxy"].map(lambda value: f"{value:.3f}")
        )
    return frame


__all__ = [
    "MTDNA_LENGTH",
    "ORIL_START",
    "MitochondrialGene",
    "MTDNA_GENES",
    "GENES_BY_KEY",
    "CANONICAL_ORDER",
    "DISPLAY_NAMES",
    "MAJOR_ARC_HEAVY_STRAND",
    "MAJOR_ARC_GENE_COLORS",
    "gene",
    "canonical_order",
    "display_names",
    "dssh_proxy",
    "arc_positions",
    "gene_metadata_frame",
]
