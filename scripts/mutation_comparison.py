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
from math import ceil

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import seaborn as sns


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

DEFAULT_GENE_LABELS = {
    "CO1": "COX1",
    "CO3": "COX3",
    "Cytb": "CytB",
    "ND2": "ND2",
    "ND6": "ND6",
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
    ax: Axes | None = None,
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
) -> SpectrumPlotResult:
    """Show a paired distribution for one mutation across 2--5 genes."""

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

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6.2, 1.45 * len(info.genes)), 5.7))
    else:
        fig = ax.figure
    x = np.arange(len(info.genes), dtype=float)

    if show_species_lines:
        for species_values in values:
            ax.plot(x, species_values, color="#718096", alpha=0.10, linewidth=0.7)
    boxplot = ax.boxplot(
        [values[:, index] for index in range(len(info.genes))],
        positions=x,
        widths=0.46,
        patch_artist=True,
        showfliers=False,
        tick_labels=[labels.get(gene, gene) for gene in info.genes],
    )
    for patch, gene in zip(boxplot["boxes"], info.genes):
        patch.set_facecolor(colors[gene])
        patch.set_alpha(0.30)
    for median in boxplot["medians"]:
        median.set_color("#1A202C")

    ordered_summary = mutation_summary.set_index(gene_col).loc[list(info.genes)]
    estimates = ordered_summary["estimate"].to_numpy(dtype=float)
    errors = np.vstack(
        [
            estimates - ordered_summary["ci_low"].to_numpy(dtype=float),
            ordered_summary["ci_high"].to_numpy(dtype=float) - estimates,
        ]
    )
    ax.errorbar(
        x,
        estimates,
        yerr=errors,
        fmt="D",
        color="#111827",
        ecolor="#111827",
        capsize=4,
        markersize=5,
        linewidth=1.4,
        label=f"Mean ({int(confidence * 100)}% bootstrap CI)",
    )
    ax.set(
        xlabel="Gene",
        ylabel="Normalized spectrum weight",
        title=f"Matched {mutation} comparison (N = {info.n_species} species)",
    )
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    return SpectrumPlotResult(
        fig, ax, mutation_summary, info.genes, info.n_species
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
        ax.set_ylabel("Mean normalized spectrum weight")
        ax.legend(frameon=False, fontsize=9)
        sns.despine(ax=ax)

    for unused_axis in flat_axes[n_panels:]:
        unused_axis.set_visible(False)
    fig.suptitle(
        f"Mutation spectrum gradient across matched species "
        f"(N = {ordered_info.n_species})",
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


def _normalise_genes(genes: Sequence[str]) -> tuple[str, ...]:
    ordered_genes = tuple(genes)
    if not 2 <= len(ordered_genes) <= 5:
        raise ValueError("genes must contain between 2 and 5 entries")
    if len(set(ordered_genes)) != len(ordered_genes):
        raise ValueError("genes must contain unique entries")
    return ordered_genes


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
        return dict(zip(genes, sns.color_palette("colorblind", len(genes))))
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
    "DEFAULT_GENE_LABELS",
    "MatchedSpectrumInfo",
    "SBS12_ORDER",
    "SpectrumPlotResult",
    "complement_substitution",
    "orient_substitutions_to_heavy_strand",
    "plot_matched_spectra",
    "plot_mutation_comparison",
    "plot_tsss_gradient",
    "select_common_species",
    "summarize_matched_spectra",
    "summarize_tsss_slopes",
    "validate_matched_spectra",
]
