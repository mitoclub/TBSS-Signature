"""Reusable plots for matched 12-component mutation spectra.

The plotting functions in this module deliberately do not calculate a species
intersection.  Their input must already contain exactly the same species for
every selected gene.  Keeping cohort construction separate makes it difficult
to compare different, silently changing species sets in one figure.

All uncertainty intervals resample whole species.  The same bootstrap draw is
therefore used for every gene and mutation component, preserving the matched
structure of the data.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations
from math import ceil
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, LogLocator
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon

_REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from mtdna import DISPLAY_NAMES, canonical_order


SBS12_ORDER = (
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
)

# Display names come from the canonical mitochondrial table so that one gene
# is never labelled differently in two stages.
DEFAULT_GENE_LABELS = dict(DISPLAY_NAMES)

DEFAULT_GENE_COLORS = {
    "CO1": "#0072B2",
    "CO3": "#E69F00",
    "Cytb": "#009E73",
    "ND2": "#CC79A7",
    "ND6": "#D55E00",
}

_OKABE_ITO = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#D55E00",
    "#56B4E9",
    "#F0E442",
    "#000000",
)
_PLOT_TEXT = "#20262E"
_PAIRED_LINE = "#667085"

DEFAULT_MUTATION_RATIOS = {
    "C>T/G>A": ("C>T", "G>A"),
    "A>G/T>C": ("A>G", "T>C"),
}

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class MatchedSpectrumInfo:
    """Validated dimensions of a matched spectrum table."""

    genes: tuple[str, ...]
    species: tuple[object, ...]
    mutations: tuple[str, ...]

    @property
    def n_species(self) -> int:
        return len(self.species)


@dataclass(frozen=True)
class SpectrumPlotResult:
    """Figure plus the numerical values used to draw it."""

    fig: Figure
    axes: Axes | np.ndarray
    summary: pd.DataFrame
    genes: tuple[str, ...]
    n_species: int
    pairwise_tests: pd.DataFrame | None = None


def complement_substitution(mutation: str) -> str:
    """Complement both alleles of a context-free substitution label."""

    parts = str(mutation).split(">")
    if len(parts) != 2 or any(
        len(base) != 1 or base not in "ACGT" for base in parts
    ):
        raise ValueError(f"Invalid single-base substitution label: {mutation!r}")
    reference, alternate = parts
    return (
        f"{reference.translate(_COMPLEMENT)}>"
        f"{alternate.translate(_COMPLEMENT)}"
    )


def orient_substitutions_to_heavy_strand(
    data: pd.DataFrame,
    *,
    mutation_col: str = "Mut",
    keep_original: bool = True,
    original_col: str = "Mut_raw",
) -> pd.DataFrame:
    """Return a copy with mutation labels complemented to heavy-strand notation.

    Use this helper only once and only when the source table is gene/coding-strand
    oriented.  Complementing an already heavy-strand-oriented table would undo
    the desired orientation.
    """

    if mutation_col not in data.columns:
        raise ValueError(f"Missing mutation column: {mutation_col!r}")
    result = data.copy()
    if keep_original:
        if original_col in result.columns and original_col != mutation_col:
            raise ValueError(
                f"Cannot preserve source labels: {original_col!r} already exists"
            )
        result[original_col] = result[mutation_col]
    result[mutation_col] = result[mutation_col].map(complement_substitution)
    return result


def select_common_species(
    data: pd.DataFrame,
    genes: Sequence[str],
    *,
    species_col: str = "Species",
    gene_col: str = "Gene",
) -> tuple[pd.DataFrame, tuple[object, ...]]:
    """Explicitly restrict a table to the all-gene species intersection.

    This is a preparation helper, not part of plotting.  It returns only the
    requested genes plus the common-species identifiers so the cohort can be
    saved or reported independently.
    """

    ordered_genes = _normalise_genes(genes)
    missing_columns = {species_col, gene_col}.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    observed_genes = set(data[gene_col].dropna())
    missing_genes = [gene for gene in ordered_genes if gene not in observed_genes]
    if missing_genes:
        raise ValueError(f"Genes absent from the input table: {missing_genes}")

    species_sets = [
        set(data.loc[data[gene_col].eq(gene), species_col].dropna())
        for gene in ordered_genes
    ]
    common_species = set.intersection(*species_sets)
    if not common_species:
        raise ValueError(f"No species are shared by all genes: {ordered_genes}")
    ordered_species = tuple(sorted(common_species, key=str))

    matched = data.loc[
        data[gene_col].isin(ordered_genes)
        & data[species_col].isin(common_species)
    ].copy()
    return matched, ordered_species


def validate_matched_spectra(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    mutation_order: Sequence[str] = SBS12_ORDER,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
    atol: float = 1e-8,
) -> MatchedSpectrumInfo:
    """Validate a complete, already matched 12-component spectrum table.

    The function is intentionally strict: it errors on extra genes, different
    species sets, duplicate components, incomplete profiles, and invalid values.
    It never intersects or drops observations on the caller's behalf.
    """

    ordered_genes = _normalise_genes(genes)
    ordered_mutations = tuple(mutation_order)
    if len(ordered_mutations) != 12 or len(set(ordered_mutations)) != 12:
        raise ValueError("mutation_order must contain 12 unique substitutions")

    required_columns = {species_col, gene_col, mutation_col, value_col}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")
    if data.empty:
        raise ValueError("The matched spectrum table is empty")
    if data[[species_col, gene_col, mutation_col, value_col]].isna().any().any():
        raise ValueError("Required identifiers and spectrum values must be complete")

    observed_genes = set(data[gene_col])
    requested_genes = set(ordered_genes)
    if observed_genes != requested_genes:
        missing = sorted(requested_genes - observed_genes, key=str)
        extra = sorted(observed_genes - requested_genes, key=str)
        raise ValueError(
            "Observed genes must exactly equal the requested genes; "
            f"missing={missing}, extra={extra}"
        )

    duplicate_mask = data.duplicated([species_col, gene_col, mutation_col])
    if duplicate_mask.any():
        example = data.loc[
            duplicate_mask, [species_col, gene_col, mutation_col]
        ].head(3)
        raise ValueError(
            "Duplicate Species/Gene/Mut rows were found; examples: "
            f"{example.to_dict('records')}"
        )

    numeric_values = pd.to_numeric(data[value_col], errors="coerce")
    if numeric_values.isna().any() or not np.isfinite(numeric_values).all():
        raise ValueError(f"{value_col} values must be numeric, complete, and finite")
    if (numeric_values < 0).any():
        raise ValueError(f"{value_col} values must be non-negative")

    expected_mutations = frozenset(ordered_mutations)
    profile_mutations = data.groupby(
        [species_col, gene_col], observed=True
    )[mutation_col].agg(frozenset)
    invalid_profiles = profile_mutations[profile_mutations.ne(expected_mutations)]
    if not invalid_profiles.empty:
        example_keys = list(invalid_profiles.index[:3])
        raise ValueError(
            f"{len(invalid_profiles)} profiles do not contain the exact 12 "
            f"substitutions; examples: {example_keys}"
        )

    species_by_gene = {
        gene: set(data.loc[data[gene_col].eq(gene), species_col])
        for gene in ordered_genes
    }
    reference_gene = ordered_genes[0]
    reference_species = species_by_gene[reference_gene]
    mismatch_messages = []
    for gene in ordered_genes[1:]:
        missing = reference_species - species_by_gene[gene]
        extra = species_by_gene[gene] - reference_species
        if missing or extra:
            mismatch_messages.append(
                f"{gene}: missing {len(missing)}, extra {len(extra)} "
                f"(examples {sorted(missing | extra, key=str)[:3]})"
            )
    if mismatch_messages:
        raise ValueError(
            "Species sets differ across genes relative to "
            f"{reference_gene}: " + "; ".join(mismatch_messages)
        )
    if not reference_species:
        raise ValueError("The common species set is empty")

    if require_sum_to_one:
        checked = data.assign(_numeric_spectrum_value=numeric_values)
        profile_sums = checked.groupby(
            [species_col, gene_col], observed=True
        )["_numeric_spectrum_value"].sum()
        invalid_sums = ~np.isclose(
            profile_sums.to_numpy(dtype=float), 1.0, rtol=0, atol=atol
        )
        if invalid_sums.any():
            examples = profile_sums.iloc[np.flatnonzero(invalid_sums)[:3]]
            raise ValueError(
                f"{value_col} must sum to one within every species/gene "
                f"profile; examples: {examples.to_dict()}"
            )

    return MatchedSpectrumInfo(
        genes=ordered_genes,
        species=tuple(sorted(reference_species, key=str)),
        mutations=ordered_mutations,
    )


def summarize_matched_spectra(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    mutation_order: Sequence[str] = SBS12_ORDER,
    estimator: str = "mean",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
) -> pd.DataFrame:
    """Summarise matched spectra with species-cluster bootstrap intervals."""

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    if estimator not in {"mean", "median"}:
        raise ValueError("estimator must be 'mean' or 'median'")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    if n_boot < 1:
        raise ValueError("n_boot must be at least one")

    values = _spectrum_array(
        data,
        info,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
    )
    estimate_function = np.mean if estimator == "mean" else np.median
    estimates = estimate_function(values, axis=0)
    bootstrap = _bootstrap_estimates(
        values,
        estimate_function=estimate_function,
        n_boot=n_boot,
        random_state=random_state,
    )
    alpha = (1 - confidence) / 2
    lows, highs = np.quantile(bootstrap, [alpha, 1 - alpha], axis=0)

    records = []
    for gene_index, gene in enumerate(info.genes):
        for mutation_index, mutation in enumerate(info.mutations):
            records.append(
                {
                    gene_col: gene,
                    mutation_col: mutation,
                    "n_species": info.n_species,
                    "estimate": float(estimates[gene_index, mutation_index]),
                    "ci_low": float(lows[gene_index, mutation_index]),
                    "ci_high": float(highs[gene_index, mutation_index]),
                    "estimator": estimator,
                    "confidence": confidence,
                }
            )
    return pd.DataFrame.from_records(records)


def plot_matched_spectra(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    mutation_order: Sequence[str] = SBS12_ORDER,
    gene_labels: Mapping[str, str] | None = None,
    kind: str = "dot",
    estimator: str = "mean",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    palette: Mapping[str, object] | Sequence[object] | None = None,
    orientation_label: str = "heavy-strand orientation",
    ax: Axes | None = None,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
) -> SpectrumPlotResult:
    """Plot estimates for 2--5 genes from an already matched DataFrame.

    ``kind='dot'`` shows bootstrap uncertainty and is the recommended default.
    ``kind='heatmap'`` is a compact overview for four or five genes; its values
    are estimates only, while the returned summary retains the intervals.
    """

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    if kind not in {"dot", "heatmap"}:
        raise ValueError("kind must be 'dot' or 'heatmap'")

    summary = summarize_matched_spectra(
        data,
        genes=info.genes,
        mutation_order=info.mutations,
        estimator=estimator,
        confidence=confidence,
        n_boot=n_boot,
        random_state=random_state,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    labels = {**DEFAULT_GENE_LABELS, **(gene_labels or {})}
    colors = _resolve_palette(info.genes, palette)

    if ax is None:
        figure_size = (14, 5.7) if kind == "dot" else (14, 2.0 + 0.55 * len(info.genes))
        fig, ax = plt.subplots(figsize=figure_size)
    else:
        fig = ax.figure

    if kind == "dot":
        x = np.arange(len(info.mutations), dtype=float)
        offsets = np.linspace(-0.28, 0.28, len(info.genes))
        for offset, gene in zip(offsets, info.genes):
            gene_summary = (
                summary.loc[summary[gene_col].eq(gene)]
                .set_index(mutation_col)
                .loc[list(info.mutations)]
            )
            estimates = gene_summary["estimate"].to_numpy(dtype=float)
            errors = np.vstack(
                [
                    estimates - gene_summary["ci_low"].to_numpy(dtype=float),
                    gene_summary["ci_high"].to_numpy(dtype=float) - estimates,
                ]
            )
            ax.errorbar(
                x + offset,
                estimates,
                yerr=errors,
                fmt="o",
                markersize=5.2,
                capsize=2.2,
                elinewidth=1.1,
                color=colors[gene],
                label=labels.get(gene, gene),
            )
        ax.set_xticks(x, info.mutations)
        ax.set_xlabel(f"Substitution ({orientation_label})")
        ax.set_ylabel(f"{estimator.capitalize()} normalized spectrum weight")
        ax.legend(title="Gene", frameon=False, ncols=min(3, len(info.genes)))
        sns.despine(ax=ax)
    else:
        matrix = (
            summary.pivot(index=gene_col, columns=mutation_col, values="estimate")
            .reindex(index=info.genes, columns=info.mutations)
        )
        matrix.index = [labels.get(gene, gene) for gene in info.genes]
        sns.heatmap(
            matrix,
            cmap="mako",
            annot=True,
            fmt=".3f",
            linewidths=0.35,
            cbar_kws={"label": f"{estimator.capitalize()} spectrum weight"},
            ax=ax,
        )
        ax.set(xlabel=f"Substitution ({orientation_label})", ylabel="Gene")

    ax.set_title(
        f"Matched 12-component mutation spectra (N = {info.n_species} species)"
    )
    fig.tight_layout()
    return SpectrumPlotResult(fig, ax, summary, info.genes, info.n_species)


def plot_mutation_comparison(
    data: pd.DataFrame,
    mutation: str,
    *,
    genes: Sequence[str],
    mutation_order: Sequence[str] = SBS12_ORDER,
    gene_labels: Mapping[str, str] | None = None,
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    palette: Mapping[str, object] | Sequence[object] | None = None,
    show_species_lines: bool = True,
    pairwise_alternative: str = "two-sided",
    ax: Axes | None = None,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
) -> SpectrumPlotResult:
    """Show a paired distribution for one mutation across 2--5 genes.

    ``pairwise_alternative`` describes Gene2 relative to Gene1 for every pair:
    ``"greater"`` tests Gene2 > Gene1, ``"less"`` tests Gene2 < Gene1, and
    ``"two-sided"`` retains the compatibility default.  The observations are
    paired by species, so the test is Wilcoxon's signed-rank test rather than
    the unpaired Mann--Whitney U test.
    """

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
    )
    if mutation not in info.mutations:
        raise ValueError(f"Unknown mutation {mutation!r}; expected {info.mutations}")
    summary = summarize_matched_spectra(
        data,
        genes=info.genes,
        mutation_order=info.mutations,
        confidence=confidence,
        n_boot=n_boot,
        random_state=random_state,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
    )
    mutation_summary = summary.loc[summary[mutation_col].eq(mutation)].copy()
    labels = {**DEFAULT_GENE_LABELS, **(gene_labels or {})}
    colors = _resolve_palette(info.genes, palette)
    values = _mutation_matrix(
        data,
        mutation,
        info,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
    )
    pairwise_tests = _paired_wilcoxon_tests(
        values,
        info.genes,
        comparison_label=mutation,
        comparison_col=mutation_col,
        alternative=pairwise_alternative,
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(3.5, 1.55 * len(info.genes)), 3.45))
    else:
        fig = ax.figure
    _draw_paired_boxplot(
        ax,
        values,
        info.genes,
        labels,
        colors,
        show_species_lines=show_species_lines,
    )
    _style_paired_axis(ax)
    ax.set_ylabel("Normalized spectrum weight")
    ax.set_title(
        f"{mutation}  ·  N = {info.n_species} paired species",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
        pad=8,
    )
    _annotate_pairwise_tests(
        ax,
        pairwise_tests,
        info.genes,
        gene_labels=labels,
        y_position=0.90,
        stacked_single=True,
    )
    fig.tight_layout(pad=0.7)
    return SpectrumPlotResult(
        fig,
        ax,
        mutation_summary,
        info.genes,
        info.n_species,
        pairwise_tests,
    )


def plot_tsss_gradient(
    data: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    *,
    genes: Sequence[str],
    mutations: Sequence[str] = ("C>T",),
    tsss_col: str = "dssh_proxy",
    label_col: str = "display_gene",
    proxy_label: str = "Relative single-stranded duration proxy (DssH)",
    mutation_order: Sequence[str] = SBS12_ORDER,
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    show_species_lines: bool = True,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
    y_label: str = "Mean normalized spectrum weight",
    figure_title: str = "Mutation spectrum gradient across matched species",
) -> SpectrumPlotResult:
    """Plot selected mutation weights against an explicit numeric TSSS proxy.

    The function does not infer biology from row or gene order.  ``gene_metadata``
    must provide one numeric proxy value per gene; points are ordered by it.
    """

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    selected_mutations = tuple(mutations)
    if not selected_mutations or len(set(selected_mutations)) != len(selected_mutations):
        raise ValueError("mutations must contain at least one unique label")
    unknown_mutations = set(selected_mutations).difference(info.mutations)
    if unknown_mutations:
        raise ValueError(f"Unknown mutations: {sorted(unknown_mutations)}")

    metadata = _validate_gene_metadata(
        gene_metadata,
        info.genes,
        gene_col=gene_col,
        tsss_col=tsss_col,
        label_col=label_col,
    )
    ordered_genes = tuple(metadata[gene_col])
    ordered_info = validate_matched_spectra(
        data,
        genes=ordered_genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    summary = summarize_matched_spectra(
        data,
        genes=ordered_genes,
        mutation_order=ordered_info.mutations,
        confidence=confidence,
        n_boot=n_boot,
        random_state=random_state,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    ).merge(metadata, on=gene_col, validate="many_to_one")

    n_panels = len(selected_mutations)
    n_columns = min(2, n_panels)
    n_rows = ceil(n_panels / n_columns)
    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(6.3 * n_columns, 4.4 * n_rows),
        sharex=True,
        squeeze=False,
    )
    flat_axes = axes.ravel()
    x = metadata[tsss_col].to_numpy(dtype=float)
    display_labels = metadata[label_col].astype(str).to_list()
    colors = dict(
        zip(selected_mutations, sns.color_palette("colorblind", n_panels))
    )

    for panel_index, mutation in enumerate(selected_mutations):
        ax = flat_axes[panel_index]
        values = _mutation_matrix(
            data,
            mutation,
            ordered_info,
            species_col=species_col,
            gene_col=gene_col,
            mutation_col=mutation_col,
            value_col=value_col,
        )
        if show_species_lines:
            for species_values in values:
                ax.plot(x, species_values, color="#64748B", alpha=0.09, linewidth=0.7)

        mutation_summary = (
            summary.loc[summary[mutation_col].eq(mutation)]
            .set_index(gene_col)
            .loc[list(ordered_genes)]
        )
        estimates = mutation_summary["estimate"].to_numpy(dtype=float)
        errors = np.vstack(
            [
                estimates - mutation_summary["ci_low"].to_numpy(dtype=float),
                mutation_summary["ci_high"].to_numpy(dtype=float) - estimates,
            ]
        )
        ax.errorbar(
            x,
            estimates,
            yerr=errors,
            color=colors[mutation],
            marker="o",
            markersize=6,
            linewidth=2,
            capsize=4,
            label=f"Mean ({int(confidence * 100)}% bootstrap CI)",
        )
        ax.set_title(mutation)
        ax.set_xticks(x, display_labels)
        ax.set_xlabel(proxy_label)
        ax.set_ylabel(y_label)
        ax.legend(frameon=False, fontsize=9)
        sns.despine(ax=ax)

    for unused_axis in flat_axes[n_panels:]:
        unused_axis.set_visible(False)
    fig.suptitle(
        f"{figure_title} (N = {ordered_info.n_species})",
        y=1.01,
    )
    fig.tight_layout()
    return SpectrumPlotResult(
        fig,
        axes if n_panels > 1 else flat_axes[0],
        summary.loc[summary[mutation_col].isin(selected_mutations)].copy(),
        ordered_genes,
        ordered_info.n_species,
    )


def summarize_tsss_slopes(
    data: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    *,
    genes: Sequence[str],
    mutations: Sequence[str] = ("C>T",),
    tsss_col: str = "dssh_proxy",
    label_col: str = "display_gene",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    mutation_order: Sequence[str] = SBS12_ORDER,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
) -> pd.DataFrame:
    """Summarise within-species linear slopes along a supplied TSSS proxy.

    These slopes are descriptive.  With only a few genes they do not isolate a
    causal TSSS effect from gene identity, sequence composition, or phylogeny.
    """

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    metadata = _validate_gene_metadata(
        gene_metadata,
        info.genes,
        gene_col=gene_col,
        tsss_col=tsss_col,
        label_col=label_col,
    )
    ordered_genes = tuple(metadata[gene_col])
    ordered_info = validate_matched_spectra(
        data,
        genes=ordered_genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    selected_mutations = tuple(mutations)
    unknown_mutations = set(selected_mutations).difference(ordered_info.mutations)
    if unknown_mutations:
        raise ValueError(f"Unknown mutations: {sorted(unknown_mutations)}")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    if n_boot < 1:
        raise ValueError("n_boot must be at least one")

    x = metadata[tsss_col].to_numpy(dtype=float)
    centered_x = x - x.mean()
    denominator = float(np.sum(centered_x**2))
    if np.isclose(denominator, 0):
        raise ValueError(f"{tsss_col} must vary across genes")
    rng = _as_rng(random_state)
    alpha = (1 - confidence) / 2
    records = []

    for mutation in selected_mutations:
        values = _mutation_matrix(
            data,
            mutation,
            ordered_info,
            species_col=species_col,
            gene_col=gene_col,
            mutation_col=mutation_col,
            value_col=value_col,
        )
        slopes = (values @ centered_x) / denominator
        draws = rng.integers(0, len(slopes), size=(n_boot, len(slopes)))
        bootstrap_means = slopes[draws].mean(axis=1)
        low, high = np.quantile(bootstrap_means, [alpha, 1 - alpha])
        records.append(
            {
                mutation_col: mutation,
                "n_species": ordered_info.n_species,
                "mean_slope_per_proxy_unit": float(slopes.mean()),
                "ci_low": float(low),
                "ci_high": float(high),
                "median_species_slope": float(np.median(slopes)),
                "fraction_species_positive": float(np.mean(slopes > 0)),
                "confidence": confidence,
            }
        )
    return pd.DataFrame.from_records(records)


def calculate_mutation_ratios(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    ratios: Mapping[str, tuple[str, str]] = DEFAULT_MUTATION_RATIOS,
    mutation_order: Sequence[str] = SBS12_ORDER,
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
    zero_denominator: str = "raise",
    species_col: str = "Species",
    class_col: str = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
) -> pd.DataFrame:
    """Calculate mutation ratios within every matched species/gene profile.

    Ratios are calculated before any cross-species summary.  The default policy
    raises on a zero denominator instead of silently adding a pseudocount.
    ``zero_denominator='drop_species'`` removes an affected species from every
    gene and ratio so that the matched design is retained.
    """

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    ratio_definitions = _normalise_ratio_definitions(ratios, info.mutations)
    if zero_denominator not in {"raise", "drop_species"}:
        raise ValueError(
            "zero_denominator must be either 'raise' or 'drop_species'"
        )

    index_columns = [species_col, gene_col]
    if class_col in data.columns:
        index_columns.insert(1, class_col)
    values = (
        data.pivot(
            index=index_columns,
            columns=mutation_col,
            values=value_col,
        )
        .reindex(columns=list(info.mutations))
        .sort_index()
    )

    invalid_species: set[object] = set()
    for _, (_, denominator_mutation) in ratio_definitions.items():
        invalid = values[denominator_mutation].le(0)
        if invalid.any():
            invalid_species.update(
                values.index.get_level_values(species_col)[invalid]
            )

    if invalid_species and zero_denominator == "raise":
        examples = sorted(invalid_species, key=str)[:5]
        raise ValueError(
            "Mutation-ratio denominators must be positive; affected species "
            f"include {examples}"
        )
    if invalid_species:
        keep = ~values.index.get_level_values(species_col).isin(invalid_species)
        values = values.loc[keep]
        if values.empty:
            raise ValueError("No matched species remain after denominator filtering")

    records: list[pd.DataFrame] = []
    for ratio_label, (numerator_mutation, denominator_mutation) in (
        ratio_definitions.items()
    ):
        ratio_table = values[[numerator_mutation, denominator_mutation]].reset_index()
        ratio_table = ratio_table.rename(
            columns={
                numerator_mutation: "NumeratorValue",
                denominator_mutation: "DenominatorValue",
            }
        )
        ratio_table["Ratio"] = ratio_label
        ratio_table["NumeratorMut"] = numerator_mutation
        ratio_table["DenominatorMut"] = denominator_mutation
        ratio_table["RatioValue"] = (
            ratio_table["NumeratorValue"] / ratio_table["DenominatorValue"]
        )
        records.append(ratio_table)

    result = pd.concat(records, ignore_index=True)
    if not np.isfinite(result["RatioValue"]).all():
        raise ValueError("Mutation ratios must be finite")
    output_columns = [species_col]
    if class_col in result.columns:
        output_columns.append(class_col)
    output_columns.extend(
        [
            gene_col,
            "Ratio",
            "NumeratorMut",
            "DenominatorMut",
            "NumeratorValue",
            "DenominatorValue",
            "RatioValue",
        ]
    )
    return result.loc[:, output_columns].sort_values(
        ["Ratio", species_col, gene_col]
    ).reset_index(drop=True)


def plot_mutation_ratio_comparison(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    ratios: Mapping[str, tuple[str, str]] = DEFAULT_MUTATION_RATIOS,
    mutation_order: Sequence[str] = SBS12_ORDER,
    gene_labels: Mapping[str, str] | None = None,
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
    zero_denominator: str = "raise",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    palette: Mapping[str, object] | Sequence[object] | None = None,
    show_species_lines: bool = True,
    pairwise_alternative: str = "two-sided",
    log_y: bool = True,
    species_col: str = "Species",
    class_col: str = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
) -> SpectrumPlotResult:
    """Show paired distributions of within-profile mutation ratios by gene.

    Each ratio is calculated separately for every species/gene profile before
    any cross-species summary. Boxplots and optional species trajectories
    therefore retain the matched-species design. ``pairwise_alternative``
    describes Gene2 relative to Gene1; see :func:`plot_mutation_comparison`.
    The log scale is the default because within-profile ratios are non-negative
    and characteristically right-skewed. If an observed numerator is exactly
    zero, a narrow linear segment at zero is combined with the logarithmic
    positive range; no pseudocount is introduced.
    """

    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    if n_boot < 1:
        raise ValueError("n_boot must be at least one")

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    ratio_definitions = _normalise_ratio_definitions(ratios, info.mutations)
    ratio_values = calculate_mutation_ratios(
        data,
        genes=info.genes,
        ratios=ratio_definitions,
        mutation_order=mutation_order,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
        zero_denominator=zero_denominator,
        species_col=species_col,
        class_col=class_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
    )
    ratio_labels = tuple(ratio_definitions)
    species = tuple(sorted(ratio_values[species_col].unique(), key=str))
    values = _ratio_value_array(
        ratio_values,
        species,
        info.genes,
        ratio_labels,
        species_col=species_col,
        gene_col=gene_col,
    )

    estimates = values.mean(axis=0)
    medians = np.median(values, axis=0)
    bootstrap = _bootstrap_estimates(
        values,
        estimate_function=np.mean,
        n_boot=n_boot,
        random_state=random_state,
    )
    alpha = (1 - confidence) / 2
    lows, highs = np.quantile(bootstrap, [alpha, 1 - alpha], axis=0)

    summary_records = []
    for gene_index, gene in enumerate(info.genes):
        for ratio_index, ratio_label in enumerate(ratio_labels):
            numerator, denominator = ratio_definitions[ratio_label]
            summary_records.append(
                {
                    gene_col: gene,
                    "Ratio": ratio_label,
                    "NumeratorMut": numerator,
                    "DenominatorMut": denominator,
                    "n_species": len(species),
                    "estimate": float(estimates[gene_index, ratio_index]),
                    "ci_low": float(lows[gene_index, ratio_index]),
                    "ci_high": float(highs[gene_index, ratio_index]),
                    "median": float(medians[gene_index, ratio_index]),
                    "confidence": confidence,
                }
            )
    summary = pd.DataFrame.from_records(summary_records)

    pairwise_test_tables = []
    for ratio_index, ratio_label in enumerate(ratio_labels):
        pairwise_test_tables.append(
            _paired_wilcoxon_tests(
                values[:, :, ratio_index],
                info.genes,
                comparison_label=ratio_label,
                comparison_col="Ratio",
                alternative=pairwise_alternative,
            )
        )
    pairwise_tests = pd.concat(pairwise_test_tables, ignore_index=True)

    labels = {**DEFAULT_GENE_LABELS, **(gene_labels or {})}
    colors = _resolve_palette(info.genes, palette)
    n_panels = len(ratio_labels)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(3.55 * n_panels + 0.15, 3.45),
        squeeze=False,
    )
    flat_axes = axes.ravel()

    for ratio_index, ratio_label in enumerate(ratio_labels):
        ax = flat_axes[ratio_index]
        ratio_matrix = values[:, :, ratio_index]
        if log_y and np.any(ratio_matrix < 0):
            raise ValueError(
                f"Log-scaled ratio panel {ratio_label!r} contains negative values"
            )
        _draw_paired_boxplot(
            ax,
            ratio_matrix,
            info.genes,
            labels,
            colors,
            show_species_lines=show_species_lines,
        )
        _style_paired_axis(ax)
        if log_y:
            if np.any(ratio_matrix == 0):
                positive_minimum = float(ratio_matrix[ratio_matrix > 0].min())
                ax.set_yscale(
                    "symlog",
                    base=10,
                    linthresh=positive_minimum / 2,
                    linscale=0.45,
                )
                ax.set_ylim(bottom=0)
                zero_panel_ticks = [
                    tick
                    for tick in (0, 0.1, 1, 10, 100, 1000)
                    if tick == 0 or tick <= float(ratio_matrix.max()) * 1.25
                ]
                ax.set_yticks(zero_panel_ticks)
            else:
                ax.set_yscale("log")
                ax.yaxis.set_major_locator(
                    LogLocator(base=10, subs=(1.0, 2.0, 5.0))
                )
            ax.yaxis.set_major_formatter(FuncFormatter(_plain_log_tick))
        ax.axhline(
            1,
            color="#9AA1AA",
            linewidth=0.8,
            linestyle=(0, (3, 3)),
            zorder=0,
        )
        _annotate_pairwise_tests(
            ax,
            pairwise_tests.loc[pairwise_tests["Ratio"].eq(ratio_label)],
            info.genes,
            gene_labels=labels,
        )
        ax.set_ylabel(
            "Within-species mutation ratio" if ratio_index == 0 else ""
        )
        ax.set_title(
            f"{chr(65 + ratio_index)}   {ratio_label.replace('/', ' / ')}",
            loc="left",
            fontsize=10.5,
            fontweight="bold",
            pad=8,
        )

    fig.suptitle(
        f"Within-species mutation ratios  ·  N = {len(species)} paired species",
        x=0.065,
        y=0.985,
        ha="left",
        fontsize=11.5,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.14, top=0.82, wspace=0.27)
    return SpectrumPlotResult(
        fig,
        axes if n_panels > 1 else flat_axes[0],
        summary,
        info.genes,
        len(species),
        pairwise_tests,
    )


def plot_tsss_ratio_gradient(
    data: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    *,
    genes: Sequence[str],
    ratios: Mapping[str, tuple[str, str]] = DEFAULT_MUTATION_RATIOS,
    tsss_col: str = "dssh_proxy",
    label_col: str = "display_gene",
    proxy_label: str = "Relative single-stranded duration proxy (DssH)",
    mutation_order: Sequence[str] = SBS12_ORDER,
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
    zero_denominator: str = "raise",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    show_species_lines: bool = True,
    species_col: str = "Species",
    class_col: str = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
) -> SpectrumPlotResult:
    """Plot per-species mutation ratios against a numeric TSSS proxy."""

    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    if n_boot < 1:
        raise ValueError("n_boot must be at least one")

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    metadata = _validate_gene_metadata(
        gene_metadata,
        info.genes,
        gene_col=gene_col,
        tsss_col=tsss_col,
        label_col=label_col,
    )
    ordered_genes = tuple(metadata[gene_col])
    ratio_definitions = _normalise_ratio_definitions(ratios, info.mutations)
    ratio_values = calculate_mutation_ratios(
        data,
        genes=ordered_genes,
        ratios=ratio_definitions,
        mutation_order=mutation_order,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
        zero_denominator=zero_denominator,
        species_col=species_col,
        class_col=class_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
    )
    ratio_labels = tuple(ratio_definitions)
    species = tuple(sorted(ratio_values[species_col].unique(), key=str))
    values = _ratio_value_array(
        ratio_values,
        species,
        ordered_genes,
        ratio_labels,
        species_col=species_col,
        gene_col=gene_col,
    )
    estimates = values.mean(axis=0)
    medians = np.median(values, axis=0)
    bootstrap = _bootstrap_estimates(
        values,
        estimate_function=np.mean,
        n_boot=n_boot,
        random_state=random_state,
    )
    alpha = (1 - confidence) / 2
    lows, highs = np.quantile(bootstrap, [alpha, 1 - alpha], axis=0)

    summary_records = []
    for gene_index, gene in enumerate(ordered_genes):
        for ratio_index, ratio_label in enumerate(ratio_labels):
            numerator, denominator = ratio_definitions[ratio_label]
            summary_records.append(
                {
                    gene_col: gene,
                    "Ratio": ratio_label,
                    "NumeratorMut": numerator,
                    "DenominatorMut": denominator,
                    "n_species": len(species),
                    "estimate": float(estimates[gene_index, ratio_index]),
                    "ci_low": float(lows[gene_index, ratio_index]),
                    "ci_high": float(highs[gene_index, ratio_index]),
                    "median": float(medians[gene_index, ratio_index]),
                    "confidence": confidence,
                }
            )
    summary = pd.DataFrame.from_records(summary_records).merge(
        metadata, on=gene_col, validate="many_to_one"
    )

    n_panels = len(ratio_labels)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(6.3 * n_panels, 4.6),
        sharex=True,
        squeeze=False,
    )
    flat_axes = axes.ravel()
    x = metadata[tsss_col].to_numpy(dtype=float)
    display_labels = metadata[label_col].astype(str).to_list()
    colors = dict(zip(ratio_labels, sns.color_palette("colorblind", n_panels)))

    for ratio_index, ratio_label in enumerate(ratio_labels):
        ax = flat_axes[ratio_index]
        ratio_matrix = values[:, :, ratio_index]
        if show_species_lines:
            for species_values in ratio_matrix:
                ax.plot(x, species_values, color="#64748B", alpha=0.09, linewidth=0.7)
        ratio_summary = (
            summary.loc[summary["Ratio"].eq(ratio_label)]
            .set_index(gene_col)
            .loc[list(ordered_genes)]
        )
        ratio_estimates = ratio_summary["estimate"].to_numpy(dtype=float)
        errors = np.vstack(
            [
                ratio_estimates - ratio_summary["ci_low"].to_numpy(dtype=float),
                ratio_summary["ci_high"].to_numpy(dtype=float) - ratio_estimates,
            ]
        )
        ax.errorbar(
            x,
            ratio_estimates,
            yerr=errors,
            color=colors[ratio_label],
            marker="o",
            markersize=6,
            linewidth=2,
            capsize=4,
            label=f"Mean ({int(confidence * 100)}% bootstrap CI)",
        )
        ax.set_title(ratio_label)
        ax.set_xticks(x, display_labels)
        ax.set_xlabel(proxy_label)
        ax.set_ylabel("Mean per-species mutation ratio")
        ax.legend(frameon=False, fontsize=9)
        sns.despine(ax=ax)

    fig.suptitle(
        f"Mutation-ratio gradient across matched species (N = {len(species)})",
        y=1.02,
    )
    fig.tight_layout()
    return SpectrumPlotResult(
        fig,
        axes if n_panels > 1 else flat_axes[0],
        summary,
        ordered_genes,
        len(species),
    )


def summarize_tsss_ratio_slopes(
    data: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    *,
    genes: Sequence[str],
    ratios: Mapping[str, tuple[str, str]] = DEFAULT_MUTATION_RATIOS,
    tsss_col: str = "dssh_proxy",
    label_col: str = "display_gene",
    mutation_order: Sequence[str] = SBS12_ORDER,
    value_col: str = "MutSpec",
    require_sum_to_one: bool = True,
    zero_denominator: str = "raise",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    species_col: str = "Species",
    class_col: str = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
) -> pd.DataFrame:
    """Summarise within-species slopes of mutation ratios along a TSSS proxy."""

    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    if n_boot < 1:
        raise ValueError("n_boot must be at least one")

    info = validate_matched_spectra(
        data,
        genes=genes,
        mutation_order=mutation_order,
        species_col=species_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
    )
    metadata = _validate_gene_metadata(
        gene_metadata,
        info.genes,
        gene_col=gene_col,
        tsss_col=tsss_col,
        label_col=label_col,
    )
    ordered_genes = tuple(metadata[gene_col])
    ratio_definitions = _normalise_ratio_definitions(ratios, info.mutations)
    ratio_values = calculate_mutation_ratios(
        data,
        genes=ordered_genes,
        ratios=ratio_definitions,
        mutation_order=mutation_order,
        value_col=value_col,
        require_sum_to_one=require_sum_to_one,
        zero_denominator=zero_denominator,
        species_col=species_col,
        class_col=class_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
    )
    ratio_labels = tuple(ratio_definitions)
    species = tuple(sorted(ratio_values[species_col].unique(), key=str))
    values = _ratio_value_array(
        ratio_values,
        species,
        ordered_genes,
        ratio_labels,
        species_col=species_col,
        gene_col=gene_col,
    )

    x = metadata[tsss_col].to_numpy(dtype=float)
    centered_x = x - x.mean()
    denominator = float(np.sum(centered_x**2))
    if np.isclose(denominator, 0):
        raise ValueError(f"{tsss_col} must vary across genes")
    slopes = np.sum(values * centered_x[None, :, None], axis=1) / denominator
    rng = _as_rng(random_state)
    draws = rng.integers(0, len(species), size=(n_boot, len(species)))
    bootstrap_means = slopes[draws].mean(axis=1)
    alpha = (1 - confidence) / 2
    lows, highs = np.quantile(bootstrap_means, [alpha, 1 - alpha], axis=0)

    records = []
    for ratio_index, ratio_label in enumerate(ratio_labels):
        numerator, denominator_mutation = ratio_definitions[ratio_label]
        ratio_slopes = slopes[:, ratio_index]
        records.append(
            {
                "Ratio": ratio_label,
                "NumeratorMut": numerator,
                "DenominatorMut": denominator_mutation,
                "n_species": len(species),
                "mean_slope_per_proxy_unit": float(ratio_slopes.mean()),
                "ci_low": float(lows[ratio_index]),
                "ci_high": float(highs[ratio_index]),
                "median_species_slope": float(np.median(ratio_slopes)),
                "fraction_species_positive": float(np.mean(ratio_slopes > 0)),
                "confidence": confidence,
            }
        )
    return pd.DataFrame.from_records(records)


def _normalise_ratio_definitions(
    ratios: Mapping[str, tuple[str, str]],
    available_mutations: Sequence[str],
) -> dict[str, tuple[str, str]]:
    if not ratios:
        raise ValueError("ratios must contain at least one definition")
    available = set(available_mutations)
    normalised: dict[str, tuple[str, str]] = {}
    for label, pair in ratios.items():
        if not str(label).strip():
            raise ValueError("ratio labels must not be empty")
        if len(pair) != 2:
            raise ValueError(f"Ratio {label!r} must contain two mutations")
        numerator, denominator = pair
        missing = {numerator, denominator}.difference(available)
        if missing:
            raise ValueError(
                f"Ratio {label!r} contains unknown mutations: {sorted(missing)}"
            )
        normalised[str(label)] = (numerator, denominator)
    return normalised


def _ratio_value_array(
    ratio_values: pd.DataFrame,
    species: Sequence[object],
    genes: Sequence[str],
    ratios: Sequence[str],
    *,
    species_col: str,
    gene_col: str,
) -> np.ndarray:
    columns = pd.MultiIndex.from_product(
        [genes, ratios], names=[gene_col, "Ratio"]
    )
    matrix = (
        ratio_values.pivot(
            index=species_col,
            columns=[gene_col, "Ratio"],
            values="RatioValue",
        )
        .reindex(index=list(species), columns=columns)
    )
    if matrix.isna().any().any():
        raise ValueError("Mutation-ratio table is not complete across species and genes")
    return matrix.to_numpy(dtype=float).reshape(
        len(species), len(genes), len(ratios)
    )


def _paired_wilcoxon_tests(
    values: np.ndarray,
    genes: Sequence[str],
    *,
    comparison_label: str,
    comparison_col: str,
    alternative: str = "two-sided",
) -> pd.DataFrame:
    """Run paired Wilcoxon tests with direction defined as Gene2 vs Gene1.

    SciPy defines its alternative on ``x - y``. The public alternatives here
    deliberately follow the more readable figure interpretation instead:
    ``"greater"`` means Gene2 > Gene1, while ``"less"`` means Gene2 < Gene1.
    """

    if values.ndim != 2 or values.shape[1] != len(genes):
        raise ValueError("values must have one column per gene")
    if not np.isfinite(values).all():
        raise ValueError("Paired-test values must be finite")
    alternative = _normalise_pairwise_alternative(alternative)
    scipy_alternative = {
        "two-sided": "two-sided",
        "greater": "less",
        "less": "greater",
    }[alternative]

    records = []
    for first_index, second_index in combinations(range(len(genes)), 2):
        first_values = values[:, first_index]
        second_values = values[:, second_index]
        differences = second_values - first_values
        if np.all(differences == 0):
            statistic, p_value = 0.0, 1.0
            p_value_underflow = False
        else:
            test = wilcoxon(
                first_values,
                second_values,
                alternative=scipy_alternative,
                zero_method="wilcox",
                method="auto",
            )
            statistic = float(test.statistic)
            raw_p_value = float(test.pvalue)
            p_value_underflow = raw_p_value == 0.0
            p_value = (
                float(np.finfo(float).tiny) if p_value_underflow else raw_p_value
            )
        first_gene = genes[first_index]
        second_gene = genes[second_index]
        if alternative == "greater":
            hypothesis = f"{second_gene} > {first_gene}"
        elif alternative == "less":
            hypothesis = f"{second_gene} < {first_gene}"
        else:
            hypothesis = f"{second_gene} != {first_gene}"
        records.append(
            {
                comparison_col: comparison_label,
                "Gene1": first_gene,
                "Gene2": second_gene,
                "n_species": values.shape[0],
                "n_nonzero_pairs": int(np.count_nonzero(differences)),
                "statistic": statistic,
                "p_value": p_value,
                "p_value_underflow": p_value_underflow,
                "median_paired_difference": float(np.median(differences)),
                "difference_definition": "Gene2 - Gene1",
                "alternative": alternative,
                "alternative_hypothesis": hypothesis,
                "scipy_alternative_on_gene1_minus_gene2": scipy_alternative,
                "test": (
                    "Paired Wilcoxon signed-rank (two-sided)"
                    if alternative == "two-sided"
                    else "Paired Wilcoxon signed-rank (one-sided)"
                ),
            }
        )

    result = pd.DataFrame.from_records(records)
    result["p_value_holm"] = _holm_adjust(result["p_value"].to_numpy(dtype=float))
    result["p_value_holm_underflow"] = result["p_value_underflow"]
    return result


def _holm_adjust(p_values: np.ndarray) -> np.ndarray:
    """Return Holm-adjusted p-values in their original order."""

    order = np.argsort(p_values)
    ordered = p_values[order]
    adjusted_ordered = np.maximum.accumulate(
        ordered * np.arange(len(ordered), 0, -1)
    )
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(adjusted_ordered, 1.0)
    return adjusted


def _annotate_pairwise_tests(
    ax: Axes,
    tests: pd.DataFrame,
    genes: Sequence[str],
    *,
    gene_labels: Mapping[str, str] | None = None,
    y_position: float = 0.97,
    stacked_single: bool = False,
) -> None:
    """Write a compact two-line paired-test result in the upper-right."""

    if tests.empty:
        return
    # Keep the existing private-call signature for compatibility. Directional
    # hypotheses remain in ``tests``; figures intentionally report only the
    # test family, sidedness, and p-value.
    _ = genes, gene_labels, stacked_single
    multiple_tests = len(tests) > 1
    alternatives = set(tests["alternative"])
    sidedness = "two-sided" if alternatives == {"two-sided"} else "one-sided"
    p_values = []
    for _, row in tests.iterrows():
        p_text = _format_p_value(
            float(row["p_value_holm"]),
            upper_bound=bool(row.get("p_value_holm_underflow", False)),
        )
        p_values.append(p_text)
    p_label = "Holm p-values" if multiple_tests else "p-value"
    lines = [
        f"Paired Wilcoxon, {sidedness}",
        f"{p_label} {', '.join(p_values)}",
    ]

    for line_number, line in enumerate(lines):
        ax.text(
            0.98,
            y_position - 0.055 * line_number,
            line,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.0,
            color=_PLOT_TEXT,
            bbox={
                "boxstyle": "square,pad=0.04",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 1.0,
            },
            zorder=10,
        )


def _format_p_value(p_value: float, *, upper_bound: bool = False) -> str:
    if not np.isfinite(p_value) or not 0 <= p_value <= 1:
        raise ValueError("p-value must be finite and lie between zero and one")
    if p_value == 0:
        upper_bound = True
        p_value = float(np.finfo(float).tiny)
    operator = "<" if upper_bound else "="
    if p_value < 0.001:
        return f"{operator} {p_value:.2e}"
    return f"{operator} {p_value:.3f}"


def _normalise_pairwise_alternative(alternative: str) -> str:
    alternative = str(alternative).casefold().strip()
    allowed = {"two-sided", "greater", "less"}
    if alternative not in allowed:
        raise ValueError(
            "pairwise_alternative must be 'two-sided', 'greater', or 'less'"
        )
    return alternative


def _draw_paired_boxplot(
    ax: Axes,
    values: np.ndarray,
    genes: Sequence[str],
    gene_labels: Mapping[str, str],
    colors: Mapping[str, object],
    *,
    show_species_lines: bool,
) -> None:
    """Draw all matched trajectories behind compact, journal-style boxes."""

    x = np.arange(len(genes), dtype=float)
    if show_species_lines:
        for species_values in values:
            ax.plot(
                x,
                species_values,
                color=_PAIRED_LINE,
                alpha=0.075,
                linewidth=0.45,
                solid_capstyle="round",
                zorder=1,
            )
    boxplot = ax.boxplot(
        [values[:, index] for index in range(len(genes))],
        positions=x,
        widths=0.38,
        whis=1.5,
        patch_artist=True,
        showfliers=False,
        tick_labels=[gene_labels.get(gene, gene) for gene in genes],
        boxprops={"linewidth": 1.0},
        whiskerprops={"color": "#4B5563", "linewidth": 0.9},
        capprops={"color": "#4B5563", "linewidth": 0.9},
        medianprops={"color": _PLOT_TEXT, "linewidth": 1.5},
        zorder=3,
    )
    for patch, gene in zip(boxplot["boxes"], genes):
        patch.set_facecolor(to_rgba(colors[gene], 0.58))
        patch.set_edgecolor(colors[gene])
    ax.set_xlim(-0.52, len(genes) - 0.48)


def _style_paired_axis(ax: Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(False)
    sns.despine(ax=ax, top=True, right=True)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#31363C")
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(axis="x", length=0, pad=6, labelsize=8.5, colors=_PLOT_TEXT)
    ax.tick_params(axis="y", length=3, width=0.8, labelsize=8, colors=_PLOT_TEXT)
    ax.yaxis.label.set_color(_PLOT_TEXT)
    ax.yaxis.label.set_size(9)


def _plain_log_tick(value: float, _position: int) -> str:
    if value == 0:
        return "0"
    if value < 0:
        return ""
    if value >= 1:
        return f"{value:g}"
    return f"{value:.2g}"


def _normalise_genes(genes: Sequence[str]) -> tuple[str, ...]:
    """Validate a gene selection and force it into canonical mtDNA order.

    Callers may list genes in any order; every figure, table, and ordered test
    then uses the fixed rCRS order from :mod:`mtdna`.
    """

    ordered_genes = tuple(genes)
    if not 2 <= len(ordered_genes) <= 5:
        raise ValueError("genes must contain between 2 and 5 entries")
    if len(set(ordered_genes)) != len(ordered_genes):
        raise ValueError("genes must contain unique entries")
    return canonical_order(ordered_genes)


def _spectrum_array(
    data: pd.DataFrame,
    info: MatchedSpectrumInfo,
    *,
    species_col: str,
    gene_col: str,
    mutation_col: str,
    value_col: str,
) -> np.ndarray:
    numeric_data = data.assign(
        _numeric_spectrum_value=pd.to_numeric(data[value_col], errors="raise")
    )
    matrix = numeric_data.pivot(
        index=species_col,
        columns=[gene_col, mutation_col],
        values="_numeric_spectrum_value",
    )
    columns = pd.MultiIndex.from_product(
        [info.genes, info.mutations], names=[gene_col, mutation_col]
    )
    matrix = matrix.reindex(index=list(info.species), columns=columns)
    return matrix.to_numpy(dtype=float).reshape(
        info.n_species, len(info.genes), len(info.mutations)
    )


def _mutation_matrix(
    data: pd.DataFrame,
    mutation: str,
    info: MatchedSpectrumInfo,
    *,
    species_col: str,
    gene_col: str,
    mutation_col: str,
    value_col: str,
) -> np.ndarray:
    matrix = (
        data.loc[data[mutation_col].eq(mutation)]
        .pivot(index=species_col, columns=gene_col, values=value_col)
        .reindex(index=list(info.species), columns=list(info.genes))
    )
    return matrix.to_numpy(dtype=float)


def _bootstrap_estimates(
    values: np.ndarray,
    *,
    estimate_function,
    n_boot: int,
    random_state: int | np.random.Generator,
    chunk_size: int = 250,
) -> np.ndarray:
    rng = _as_rng(random_state)
    bootstrap = np.empty((n_boot, values.shape[1], values.shape[2]), dtype=float)
    for start in range(0, n_boot, chunk_size):
        stop = min(start + chunk_size, n_boot)
        indices = rng.integers(
            0, values.shape[0], size=(stop - start, values.shape[0])
        )
        bootstrap[start:stop] = estimate_function(values[indices], axis=1)
    return bootstrap


def _as_rng(
    random_state: int | np.random.Generator,
) -> np.random.Generator:
    if isinstance(random_state, np.random.Generator):
        return random_state
    return np.random.default_rng(random_state)


def _resolve_palette(
    genes: Sequence[str],
    palette: Mapping[str, object] | Sequence[object] | None,
) -> dict[str, object]:
    if palette is None:
        resolved: dict[str, object] = {}
        fallback_index = 0
        for gene in genes:
            if gene in DEFAULT_GENE_COLORS:
                resolved[gene] = DEFAULT_GENE_COLORS[gene]
            else:
                resolved[gene] = _OKABE_ITO[fallback_index % len(_OKABE_ITO)]
                fallback_index += 1
        return resolved
    if isinstance(palette, Mapping):
        missing = [gene for gene in genes if gene not in palette]
        if missing:
            raise ValueError(f"Palette is missing colors for: {missing}")
        return {gene: palette[gene] for gene in genes}
    colors = list(palette)
    if len(colors) < len(genes):
        raise ValueError("Palette sequence has fewer colors than genes")
    return dict(zip(genes, colors))


def _validate_gene_metadata(
    gene_metadata: pd.DataFrame,
    genes: Sequence[str],
    *,
    gene_col: str,
    tsss_col: str,
    label_col: str,
) -> pd.DataFrame:
    required = {gene_col, tsss_col, label_col}
    missing = required.difference(gene_metadata.columns)
    if missing:
        raise ValueError(f"gene_metadata is missing columns: {sorted(missing)}")
    metadata = gene_metadata.loc[
        gene_metadata[gene_col].isin(genes), [gene_col, tsss_col, label_col]
    ].copy()
    if metadata[gene_col].duplicated().any():
        raise ValueError("gene_metadata must contain one row per gene")
    if set(metadata[gene_col]) != set(genes):
        missing_genes = sorted(set(genes) - set(metadata[gene_col]), key=str)
        raise ValueError(f"gene_metadata is missing genes: {missing_genes}")
    metadata[tsss_col] = pd.to_numeric(metadata[tsss_col], errors="coerce")
    if metadata[[tsss_col, label_col]].isna().any().any():
        raise ValueError("TSSS proxy values and display labels must be complete")
    if not np.isfinite(metadata[tsss_col]).all():
        raise ValueError("TSSS proxy values must be finite")
    if metadata[tsss_col].duplicated().any():
        raise ValueError("TSSS proxy values must be unique across genes")
    return metadata.sort_values(tsss_col).reset_index(drop=True)


__all__ = [
    "DEFAULT_GENE_COLORS",
    "DEFAULT_GENE_LABELS",
    "DEFAULT_MUTATION_RATIOS",
    "MatchedSpectrumInfo",
    "SBS12_ORDER",
    "SpectrumPlotResult",
    "calculate_mutation_ratios",
    "complement_substitution",
    "orient_substitutions_to_heavy_strand",
    "plot_matched_spectra",
    "plot_mutation_ratio_comparison",
    "plot_mutation_comparison",
    "plot_tsss_gradient",
    "plot_tsss_ratio_gradient",
    "select_common_species",
    "summarize_matched_spectra",
    "summarize_tsss_slopes",
    "summarize_tsss_ratio_slopes",
    "validate_matched_spectra",
]
