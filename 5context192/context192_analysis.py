"""Reusable analysis helpers for matched 192-component mutation spectra.

The source spectra use labels such as ``A[C>T]G``: the two outside bases are
the immediate neighbours and the substitution is inside brackets.  Functions
in this module keep the complete 192-component normalization explicit.  In
particular, summing the 16 context weights for a central substitution is a
useful *marginal of the 192-component model*, but it must not be presented as
the independently derived 12-component spectrum.

All bootstrap intervals resample whole species.  For complete matched spectra,
one species draw is shared by every gene and all 192 components in a replicate,
so profile covariance and mean-spectrum normalization are preserved.
``weighting='species'`` gives every species the same weight.
``weighting='equal_class'`` first estimates a quantity within each represented
vertebrate class and then averages those class estimates equally; this is a
sensitivity analysis, not a phylogenetic correction.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import re

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
import pandas as pd


BASE_ORDER = ("A", "C", "G", "T")

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

# PyMutSpec's ``ordered_sbs192_kp`` groups strand-complementary substitutions
# next to one another.  The paired order and colors below reproduce the visual
# grammar used by the reference mtdna-192component-mutspec-chordata project.
PAIRED_SBS12_ORDER = (
    "C>A",
    "G>T",
    "C>G",
    "G>C",
    "C>T",
    "G>A",
    "T>A",
    "A>T",
    "T>C",
    "A>G",
    "T>G",
    "A>C",
)

SBS12_COLOR_MAPPING = {
    "C>A": "deepskyblue",
    "G>T": "deepskyblue",
    "C>G": "black",
    "G>C": "black",
    "C>T": "red",
    "G>A": "red",
    "T>A": "silver",
    "A>T": "silver",
    "T>C": "yellowgreen",
    "A>G": "yellowgreen",
    "T>G": "pink",
    "A>C": "pink",
}

_COSMIC_SBS6 = frozenset(("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# The second member of each pair is ordered by the reverse-complemented
# context, as in PyMutSpec.  This aligns mirrored contexts across paired
# substitution blocks instead of sorting every block independently.
PAIRED_CONTEXT192_ORDER = tuple(
    (
        f"{left}[{substitution}]{right}"
        if substitution in _COSMIC_SBS6
        else (
            f"{right.translate(_COMPLEMENT)}[{substitution}]"
            f"{left.translate(_COMPLEMENT)}"
        )
    )
    for substitution in PAIRED_SBS12_ORDER
    for left in BASE_ORDER
    for right in BASE_ORDER
)

CONTEXT192_ORDER = tuple(
    f"{left}[{substitution}]{right}"
    for substitution in SBS12_ORDER
    for left in BASE_ORDER
    for right in BASE_ORDER
)

QC_STATUS_ORDER = (
    "no_expected_opportunity",
    "opportunity_observed_zero",
    "opportunity_observed_positive",
    "no_expected_opportunity_but_observed_positive",
)

_CONTEXT_PATTERN = re.compile(
    r"^(?P<left>[ACGT])\[(?P<reference>[ACGT])>"
    r"(?P<alternate>[ACGT])\](?P<right>[ACGT])$"
)


@dataclass(frozen=True)
class ContextMutation:
    """Parsed representation of one ``L[REF>ALT]R`` component."""

    left: str
    reference: str
    alternate: str
    right: str

    @property
    def central_substitution(self) -> str:
        return f"{self.reference}>{self.alternate}"

    @property
    def flanking_context(self) -> str:
        return f"{self.left}_{self.right}"

    @property
    def reference_trinucleotide(self) -> str:
        return f"{self.left}{self.reference}{self.right}"

    def __str__(self) -> str:
        return (
            f"{self.left}[{self.reference}>{self.alternate}]"
            f"{self.right}"
        )


@dataclass(frozen=True)
class Context192ProfileInfo:
    """Dimensions and row-level opportunity QC for a matched table."""

    genes: tuple[str, ...]
    species: tuple[object, ...]
    contexts: tuple[str, ...]
    classes: tuple[object, ...]
    class_counts: tuple[tuple[object, int], ...]
    n_no_expected_opportunity: int
    n_opportunity_observed_zero: int
    n_opportunity_observed_positive: int

    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_classes(self) -> int:
        return len(self.classes)


def parse_context_mutation(label: str) -> ContextMutation:
    """Parse and validate an exact ``L[REF>ALT]R`` mutation label."""

    match = _CONTEXT_PATTERN.fullmatch(str(label))
    if match is None:
        raise ValueError(
            f"Invalid context mutation label {label!r}; expected L[REF>ALT]R "
            "with A/C/G/T bases"
        )
    parsed = ContextMutation(**match.groupdict())
    if parsed.reference == parsed.alternate:
        raise ValueError(
            f"Invalid context mutation label {label!r}: reference and "
            "alternate bases must differ"
        )
    return parsed


def reverse_complement_context(label: str) -> str:
    """Reverse-complement a complete context and its central substitution.

    ``L[REF>ALT]R`` becomes
    ``comp(R)[comp(REF)>comp(ALT)]comp(L)``.  Swapping the flanks is essential;
    complementing only the central substitution is not a valid context
    orientation change.
    """

    parsed = parse_context_mutation(label)
    return str(
        ContextMutation(
            left=parsed.right.translate(_COMPLEMENT),
            reference=parsed.reference.translate(_COMPLEMENT),
            alternate=parsed.alternate.translate(_COMPLEMENT),
            right=parsed.left.translate(_COMPLEMENT),
        )
    )


def orient_contexts_to_heavy_strand(
    data: pd.DataFrame,
    *,
    mutation_col: str = "Mut",
    keep_original: bool = True,
    original_col: str = "Mut_raw",
    orientation_col: str = "ContextOrientation",
) -> pd.DataFrame:
    """Reverse-complement coding/light-strand labels exactly once.

    The returned table carries ``ContextOrientation='heavy'``.  Calling this
    function again on that table raises instead of silently undoing the first
    transformation.  When a source already has an orientation column, its only
    accepted pre-conversion value is ``'coding'``.  An unmarked source is
    assumed to use the coding/light-strand convention of the supplied dataset.
    """

    if mutation_col not in data.columns:
        raise ValueError(f"Missing mutation column: {mutation_col!r}")
    if orientation_col in data.columns:
        orientation = set(data[orientation_col].dropna().astype(str).str.lower())
        if orientation == {"heavy"}:
            raise ValueError("Contexts are already marked as heavy-strand oriented")
        if orientation != {"coding"}:
            raise ValueError(
                f"{orientation_col!r} must contain only 'coding' before "
                f"orientation; observed={sorted(orientation)}"
            )
    if keep_original and original_col != mutation_col and original_col in data.columns:
        raise ValueError(
            f"Cannot preserve source labels: {original_col!r} already exists"
        )

    unique_labels = pd.unique(data[mutation_col])
    if pd.isna(unique_labels).any():
        raise ValueError(f"{mutation_col!r} must not contain missing labels")
    conversion = {
        label: reverse_complement_context(str(label)) for label in unique_labels
    }

    result = data.copy()
    if keep_original and original_col != mutation_col:
        result[original_col] = result[mutation_col]
    result[mutation_col] = result[mutation_col].map(conversion)
    result[orientation_col] = "heavy"
    result.attrs["context_orientation"] = "heavy"
    return result


def select_common_species_192(
    data: pd.DataFrame,
    genes: Sequence[str],
    *,
    species_col: str = "Species",
    gene_col: str = "Gene",
) -> tuple[pd.DataFrame, tuple[object, ...]]:
    """Restrict a table to species present for every requested gene.

    No taxonomic class is excluded.  This helper constructs the intersection;
    call :func:`validate_matched_context192` afterwards to require complete,
    normalized 192-component profiles.
    """

    ordered_genes = _normalise_genes(genes)
    required = {species_col, gene_col}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if data[[species_col, gene_col]].isna().any().any():
        raise ValueError("Species and gene identifiers must not be missing")

    available_genes = set(data[gene_col])
    absent = [gene for gene in ordered_genes if gene not in available_genes]
    if absent:
        raise ValueError(f"Genes absent from the input table: {absent}")
    species_sets = [
        set(data.loc[data[gene_col].eq(gene), species_col])
        for gene in ordered_genes
    ]
    common = set.intersection(*species_sets)
    if not common:
        raise ValueError(f"No species are shared by all genes: {ordered_genes}")
    species = tuple(sorted(common, key=str))
    matched = data.loc[
        data[gene_col].isin(ordered_genes) & data[species_col].isin(common)
    ].copy()
    return matched, species


def annotate_context_qc(
    data: pd.DataFrame,
    *,
    expected_col: str = "Expected",
    observed_col: str = "Observed",
    status_col: str = "OpportunityStatus",
) -> pd.DataFrame:
    """Annotate expected-opportunity and observed-event states.

    Missing ``Observed`` values in the supplied spectra encode zero inferred
    events and are therefore grouped with explicit zeros.  They remain missing
    in the original column.  ``Expected == 0`` is labelled as no expected
    opportunity under this synonymous-site derivation; it is not evidence that
    the biochemical mutation itself is impossible.
    """

    missing = {expected_col, observed_col}.difference(data.columns)
    if missing:
        raise ValueError(f"Missing QC columns: {sorted(missing)}")

    expected = pd.to_numeric(data[expected_col], errors="coerce")
    if expected.isna().any() or not np.isfinite(expected).all():
        raise ValueError(f"{expected_col} must be numeric, complete, and finite")
    if expected.lt(0).any():
        raise ValueError(f"{expected_col} must be non-negative")

    observed_raw = data[observed_col]
    observed = pd.to_numeric(observed_raw, errors="coerce")
    invalid_observed = observed_raw.notna() & observed.isna()
    if invalid_observed.any():
        raise ValueError(f"Non-missing {observed_col} values must be numeric")
    finite_observed = observed.dropna()
    if not np.isfinite(finite_observed).all() or finite_observed.lt(0).any():
        raise ValueError(
            f"Non-missing {observed_col} values must be finite and non-negative"
        )

    no_opportunity = expected.eq(0)
    observed_positive = observed.fillna(0).gt(0)
    statuses = np.select(
        [
            no_opportunity & ~observed_positive,
            ~no_opportunity & ~observed_positive,
            ~no_opportunity & observed_positive,
            no_opportunity & observed_positive,
        ],
        list(QC_STATUS_ORDER),
        default="unclassified",
    )

    result = data.copy()
    result[status_col] = pd.Categorical(
        statuses, categories=list(QC_STATUS_ORDER), ordered=True
    )
    result["HasExpectedOpportunity"] = ~no_opportunity
    result["HasObservedMutation"] = observed_positive
    return result


def validate_matched_context192(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    context_order: Sequence[str] = CONTEXT192_ORDER,
    species_col: str = "Species",
    class_col: str | None = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    value_col: str = "MutSpec",
    expected_col: str = "Expected",
    observed_col: str = "Observed",
    require_sum_to_one: bool = True,
    atol: float = 1e-8,
) -> Context192ProfileInfo:
    """Validate exact, normalized, all-gene matched 192-component profiles.

    Validation distinguishes rows with no expected synonymous opportunity from
    rows that have opportunity but no inferred observed event.  Both should
    have zero ``MutSpec``, but they have different interpretations.
    """

    ordered_genes = _normalise_genes(genes)
    contexts = _normalise_context_order(context_order)
    required = {
        species_col,
        gene_col,
        mutation_col,
        value_col,
        expected_col,
        observed_col,
    }
    if class_col is not None:
        required.add(class_col)
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if data.empty:
        raise ValueError("The matched 192-component table is empty")

    identifier_columns = [species_col, gene_col, mutation_col]
    if class_col is not None:
        identifier_columns.append(class_col)
    if data[identifier_columns].isna().any().any():
        raise ValueError("Required identifiers must not be missing")

    observed_genes = set(data[gene_col])
    requested_genes = set(ordered_genes)
    if observed_genes != requested_genes:
        missing_genes = sorted(requested_genes - observed_genes, key=str)
        extra_genes = sorted(observed_genes - requested_genes, key=str)
        raise ValueError(
            "Observed genes must exactly equal requested genes; "
            f"missing={missing_genes}, extra={extra_genes}"
        )

    duplicate = data.duplicated([species_col, gene_col, mutation_col])
    if duplicate.any():
        examples = data.loc[
            duplicate, [species_col, gene_col, mutation_col]
        ].head(3)
        raise ValueError(
            "Duplicate Species/Gene/Mut rows were found; examples: "
            f"{examples.to_dict('records')}"
        )

    values = pd.to_numeric(data[value_col], errors="coerce")
    if values.isna().any() or not np.isfinite(values).all():
        raise ValueError(f"{value_col} must be numeric, complete, and finite")
    if values.lt(0).any():
        raise ValueError(f"{value_col} must be non-negative")

    observed_contexts = set(data[mutation_col].astype(str))
    expected_contexts = set(contexts)
    if not observed_contexts.issubset(expected_contexts):
        invalid = sorted(observed_contexts - expected_contexts)[:5]
        raise ValueError(f"Unknown 192-component labels: {invalid}")
    profile_contexts = data.groupby(
        [species_col, gene_col], observed=True, sort=False
    )[mutation_col].agg(frozenset)
    expected_set = frozenset(contexts)
    invalid_profiles = profile_contexts[profile_contexts.ne(expected_set)]
    if not invalid_profiles.empty:
        examples = list(invalid_profiles.index[:3])
        raise ValueError(
            f"{len(invalid_profiles)} profiles do not contain the exact 192 "
            f"components; examples: {examples}"
        )

    species_by_gene = {
        gene: set(data.loc[data[gene_col].eq(gene), species_col])
        for gene in ordered_genes
    }
    reference_gene = ordered_genes[0]
    species = species_by_gene[reference_gene]
    mismatches = []
    for gene in ordered_genes[1:]:
        missing_species = species - species_by_gene[gene]
        extra_species = species_by_gene[gene] - species
        if missing_species or extra_species:
            mismatches.append(
                f"{gene}: missing={len(missing_species)}, "
                f"extra={len(extra_species)}"
            )
    if mismatches:
        raise ValueError(
            "Species sets differ across genes relative to "
            f"{reference_gene}: {'; '.join(mismatches)}"
        )
    if not species:
        raise ValueError("The common species set is empty")

    if require_sum_to_one:
        profile_sums = data.assign(_numeric_value=values).groupby(
            [species_col, gene_col], observed=True, sort=False
        )["_numeric_value"].sum()
        invalid_sum = ~np.isclose(
            profile_sums.to_numpy(dtype=float), 1.0, rtol=0, atol=atol
        )
        if invalid_sum.any():
            examples = profile_sums.iloc[np.flatnonzero(invalid_sum)[:3]]
            raise ValueError(
                f"{value_col} must sum to one in every species/gene "
                f"192-profile; examples: {examples.to_dict()}"
            )

    qc = annotate_context_qc(
        data, expected_col=expected_col, observed_col=observed_col
    )
    invalid_status = qc["OpportunityStatus"].eq(
        "no_expected_opportunity_but_observed_positive"
    )
    if invalid_status.any():
        raise ValueError(
            "Positive observed values were found where Expected == 0"
        )
    no_opportunity = qc["OpportunityStatus"].eq("no_expected_opportunity")
    observed_zero = qc["OpportunityStatus"].eq(
        "opportunity_observed_zero"
    )
    if values.loc[no_opportunity].gt(atol).any():
        raise ValueError(f"{value_col} must be zero where Expected == 0")
    if values.loc[observed_zero].gt(atol).any():
        raise ValueError(
            f"{value_col} must be zero where Expected > 0 but no event was observed"
        )

    if class_col is None:
        classes: tuple[object, ...] = ()
        class_counts: tuple[tuple[object, int], ...] = ()
    else:
        species_classes = data[[species_col, class_col]].drop_duplicates()
        ambiguous = species_classes.groupby(species_col, observed=True)[
            class_col
        ].nunique()
        if ambiguous.gt(1).any():
            examples = list(ambiguous[ambiguous.gt(1)].index[:3])
            raise ValueError(
                f"Species map to multiple {class_col} values; examples: {examples}"
            )
        counts = species_classes[class_col].value_counts(sort=False)
        classes = tuple(sorted(counts.index, key=str))
        class_counts = tuple((item, int(counts[item])) for item in classes)

    status_counts = qc["OpportunityStatus"].value_counts()
    return Context192ProfileInfo(
        genes=ordered_genes,
        species=tuple(sorted(species, key=str)),
        contexts=contexts,
        classes=classes,
        class_counts=class_counts,
        n_no_expected_opportunity=int(
            status_counts.get("no_expected_opportunity", 0)
        ),
        n_opportunity_observed_zero=int(
            status_counts.get("opportunity_observed_zero", 0)
        ),
        n_opportunity_observed_positive=int(
            status_counts.get("opportunity_observed_positive", 0)
        ),
    )


def add_context_columns(
    data: pd.DataFrame,
    *,
    mutation_col: str = "Mut",
) -> pd.DataFrame:
    """Return a copy with parsed central and flanking context columns."""

    if mutation_col not in data.columns:
        raise ValueError(f"Missing mutation column: {mutation_col!r}")
    labels = pd.unique(data[mutation_col])
    if pd.isna(labels).any():
        raise ValueError(f"{mutation_col!r} must not contain missing labels")
    parsed = {label: parse_context_mutation(str(label)) for label in labels}
    result = data.copy()
    result["Left"] = result[mutation_col].map(
        {label: item.left for label, item in parsed.items()}
    )
    result["Reference"] = result[mutation_col].map(
        {label: item.reference for label, item in parsed.items()}
    )
    result["Alternate"] = result[mutation_col].map(
        {label: item.alternate for label, item in parsed.items()}
    )
    result["Right"] = result[mutation_col].map(
        {label: item.right for label, item in parsed.items()}
    )
    result["CentralSubstitution"] = result[mutation_col].map(
        {label: item.central_substitution for label, item in parsed.items()}
    )
    result["FlankingContext"] = result[mutation_col].map(
        {label: item.flanking_context for label, item in parsed.items()}
    )
    result["ReferenceTrinucleotide"] = result[mutation_col].map(
        {label: item.reference_trinucleotide for label, item in parsed.items()}
    )
    return result


def calculate_within_substitution_shares(
    data: pd.DataFrame,
    *,
    value_col: str = "MutSpec",
    output_col: str = "WithinSubstitutionShare",
    total_col: str = "CentralWeight192",
    zero_total: str = "nan",
    species_col: str = "Species",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
) -> pd.DataFrame:
    """Calculate context shares conditional on each central substitution.

    ``CentralWeight192`` is the sum of 16 components *within the normalized
    192-component spectrum*.  It is deliberately not named ``MutSpec12``:
    aggregation after context-specific opportunity correction is not generally
    equal to a separately derived 12-component spectrum.

    For a species/gene/substitution whose 16 components sum to zero, the share
    is undefined.  ``zero_total='nan'`` preserves those rows as missing,
    ``'drop'`` removes them, and ``'raise'`` stops immediately.  No pseudocount
    is added.
    """

    if zero_total not in {"nan", "drop", "raise"}:
        raise ValueError("zero_total must be 'nan', 'drop', or 'raise'")
    required = {species_col, gene_col, mutation_col, value_col}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    result = add_context_columns(data, mutation_col=mutation_col)
    values = pd.to_numeric(result[value_col], errors="coerce")
    if values.isna().any() or not np.isfinite(values).all():
        raise ValueError(f"{value_col} must be numeric, complete, and finite")
    if values.lt(0).any():
        raise ValueError(f"{value_col} must be non-negative")
    result[value_col] = values
    group_columns = [species_col, gene_col, "CentralSubstitution"]
    result[total_col] = result.groupby(
        group_columns, observed=True, sort=False
    )[value_col].transform("sum")
    zero = result[total_col].eq(0)
    if zero.any() and zero_total == "raise":
        examples = result.loc[zero, group_columns].drop_duplicates().head(5)
        raise ValueError(
            "Within-substitution shares are undefined for zero-total groups; "
            f"examples: {examples.to_dict('records')}"
        )
    result[output_col] = result[value_col] / result[total_col].where(~zero)
    if zero_total == "drop":
        result = result.loc[~zero].copy()
    return result


def summarize_context_components(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    central_substitutions: Sequence[str] | None = None,
    contexts: Sequence[str] | None = None,
    value_col: str = "MutSpec",
    estimator: str = "mean",
    weighting: str = "species",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    require_complete_profiles: bool = True,
    context_order: Sequence[str] = CONTEXT192_ORDER,
    species_col: str = "Species",
    class_col: str = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    expected_col: str = "Expected",
    observed_col: str = "Observed",
) -> pd.DataFrame:
    """Summarize context components with species-bootstrap intervals.

    With ``require_complete_profiles=True`` (the default), the input must pass
    exact 192-profile validation.  One matched-species sample is then drawn per
    bootstrap replicate and reused for every selected gene and context.  This
    preserves covariance across the 192 components and across matched genes;
    for mean complete spectra, every bootstrap spectrum still sums to one.

    Set ``require_complete_profiles=False`` for a deliberate incomplete subset
    such as ``Expected > 0`` opportunity-normalized O/E rows.  Because the
    available species can then differ by component, each gene/context group is
    bootstrapped independently and its own ``n_species`` is reported.  Missing
    rows are never converted to zero.
    """

    ordered_genes = _normalise_genes(genes)
    ordered_contexts = _normalise_context_order(context_order)
    _validate_bootstrap_options(estimator, weighting, confidence, n_boot)
    if weighting == "equal_class" and class_col not in data.columns:
        raise ValueError(f"weighting='equal_class' requires {class_col!r}")

    if require_complete_profiles:
        validate_matched_context192(
            data,
            genes=ordered_genes,
            context_order=ordered_contexts,
            species_col=species_col,
            class_col=class_col,
            gene_col=gene_col,
            mutation_col=mutation_col,
            value_col=value_col,
            expected_col=expected_col,
            observed_col=observed_col,
        )
    prepared = _prepare_component_rows(
        data,
        genes=ordered_genes,
        value_col=value_col,
        allow_missing_values=not require_complete_profiles,
        species_col=species_col,
        class_col=class_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
    )
    selected_contexts = _select_contexts(
        ordered_contexts,
        central_substitutions=central_substitutions,
        contexts=contexts,
    )
    prepared = prepared.loc[prepared[mutation_col].isin(selected_contexts)]
    if prepared.empty:
        raise ValueError("No rows remain for the selected contexts")

    estimate_function = np.mean if estimator == "mean" else np.median
    rng = _as_rng(random_state)
    alpha = (1 - confidence) / 2
    records: list[dict[str, object]] = []
    if require_complete_profiles:
        (
            species,
            species_classes,
            joint_estimate,
            joint_bootstrap,
        ) = _bootstrap_matched_profile_estimate(
            prepared,
            genes=ordered_genes,
            contexts=selected_contexts,
            value_col=value_col,
            estimator=estimate_function,
            weighting=weighting,
            n_boot=n_boot,
            rng=rng,
            species_col=species_col,
            class_col=class_col,
            gene_col=gene_col,
            mutation_col=mutation_col,
        )
        joint_low = np.quantile(joint_bootstrap, alpha, axis=0)
        joint_high = np.quantile(joint_bootstrap, 1 - alpha, axis=0)
        n_classes = int(pd.Series(species_classes).nunique())
        for gene_index, gene in enumerate(ordered_genes):
            for context_index, context in enumerate(selected_contexts):
                estimate = float(joint_estimate[gene_index, context_index])
                low = float(joint_low[gene_index, context_index])
                high = float(joint_high[gene_index, context_index])
                parsed = parse_context_mutation(context)
                records.append(
                    {
                        gene_col: gene,
                        mutation_col: context,
                        "Left": parsed.left,
                        "Reference": parsed.reference,
                        "Alternate": parsed.alternate,
                        "Right": parsed.right,
                        "CentralSubstitution": parsed.central_substitution,
                        "FlankingContext": parsed.flanking_context,
                        "ReferenceTrinucleotide": parsed.reference_trinucleotide,
                        "n_species": int(len(species)),
                        "n_classes": n_classes,
                        "estimate": estimate,
                        "ci_low": low,
                        "ci_high": high,
                        "estimator": estimator,
                        "weighting": weighting,
                        "confidence": confidence,
                        "n_boot": n_boot,
                        "bootstrap_unit": "matched_species_profile",
                        "joint_profile_bootstrap": True,
                        "value_col": value_col,
                    }
                )
        return pd.DataFrame.from_records(records)

    for gene in ordered_genes:
        gene_data = prepared.loc[prepared[gene_col].eq(gene)]
        for context in selected_contexts:
            component = gene_data.loc[gene_data[mutation_col].eq(context)]
            if component.empty:
                continue
            values = component[value_col].to_numpy(dtype=float)
            classes = component[class_col].to_numpy()
            estimate, bootstrap = _bootstrap_group_estimate(
                values,
                classes,
                estimator=estimate_function,
                weighting=weighting,
                n_boot=n_boot,
                rng=rng,
            )
            low, high = np.quantile(bootstrap, [alpha, 1 - alpha])
            parsed = parse_context_mutation(context)
            records.append(
                {
                    gene_col: gene,
                    mutation_col: context,
                    "Left": parsed.left,
                    "Reference": parsed.reference,
                    "Alternate": parsed.alternate,
                    "Right": parsed.right,
                    "CentralSubstitution": parsed.central_substitution,
                    "FlankingContext": parsed.flanking_context,
                    "ReferenceTrinucleotide": parsed.reference_trinucleotide,
                    "n_species": int(component[species_col].nunique()),
                    "n_classes": int(component[class_col].nunique()),
                    "estimate": float(estimate),
                    "ci_low": float(low),
                    "ci_high": float(high),
                    "estimator": estimator,
                    "weighting": weighting,
                    "confidence": confidence,
                    "n_boot": n_boot,
                    "bootstrap_unit": "species_within_gene_context",
                    "joint_profile_bootstrap": False,
                    "value_col": value_col,
                }
            )
    return pd.DataFrame.from_records(records)


def plot_context192_bar_spectra(
    summary: pd.DataFrame,
    *,
    genes: Sequence[str],
    gene_labels: Mapping[str, str] | None = None,
    context_order: Sequence[str] = PAIRED_CONTEXT192_ORDER,
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    estimate_col: str = "estimate",
    ci_low_col: str | None = "ci_low",
    ci_high_col: str | None = "ci_high",
    n_species_col: str | None = "n_species",
    title: str = "Matched 192-component mutation spectra",
    ylabel: str = "Mean normalized 192-component weight",
    figsize: tuple[float, float] | None = None,
    ticksize: float = 5.5,
    titlesize: float = 16,
    fontname: str = "DejaVu Sans Mono",
    group_gap: float = 1.5,
) -> tuple[Figure, np.ndarray]:
    """Plot reference-style 192-component bars for one or more genes.

    This reproduces the visual grammar of PyMutSpec's ``plot_mutspec192``:
    strand-complementary substitutions are adjacent, each substitution forms
    a 16-context color block, paired blocks share a color, and mirrored blocks
    use reverse-complement-aligned context ordering.  Unlike the legacy helper,
    this function consumes an explicit summary table and can draw its stored
    species-bootstrap confidence intervals without importing PyMutSpec.

    ``summary`` must contain exactly one row for every selected gene/context.
    It is expected to come from :func:`summarize_context_components`; the
    function never reorients mutation labels or recalculates the cohort.
    """

    ordered_genes = tuple(genes)
    if not ordered_genes or len(set(ordered_genes)) != len(ordered_genes):
        raise ValueError("genes must contain at least one unique entry")
    ordered_contexts = _normalise_context_order(context_order)
    if tuple(
        parse_context_mutation(context).central_substitution
        for context in ordered_contexts[::16]
    ) != PAIRED_SBS12_ORDER:
        raise ValueError(
            "context_order must contain 16 consecutive contexts for each "
            "substitution in PAIRED_SBS12_ORDER"
        )

    required = {gene_col, mutation_col, estimate_col}
    use_intervals = ci_low_col is not None or ci_high_col is not None
    if use_intervals:
        if ci_low_col is None or ci_high_col is None:
            raise ValueError(
                "ci_low_col and ci_high_col must both be supplied or both be None"
            )
        required.update((ci_low_col, ci_high_col))
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"Missing plotting columns: {sorted(missing)}")

    plot_data = summary.loc[summary[gene_col].isin(ordered_genes)].copy()
    observed_genes = set(plot_data[gene_col])
    if observed_genes != set(ordered_genes):
        absent = sorted(set(ordered_genes) - observed_genes, key=str)
        raise ValueError(f"Summary is missing genes: {absent}")
    if plot_data.duplicated([gene_col, mutation_col]).any():
        raise ValueError("Summary contains duplicate Gene/Mut rows")

    expected_context_set = set(ordered_contexts)
    contexts_by_gene = plot_data.groupby(gene_col, observed=True)[
        mutation_col
    ].agg(set)
    invalid_genes = [
        gene
        for gene in ordered_genes
        if contexts_by_gene.get(gene, set()) != expected_context_set
    ]
    if invalid_genes:
        raise ValueError(
            "Every plotted gene must contain exactly 192 contexts; invalid "
            f"genes: {invalid_genes}"
        )

    numeric_columns = [estimate_col]
    if use_intervals:
        numeric_columns.extend((ci_low_col, ci_high_col))
    for column in numeric_columns:
        plot_data[column] = pd.to_numeric(plot_data[column], errors="coerce")
        if plot_data[column].isna().any() or not np.isfinite(
            plot_data[column]
        ).all():
            raise ValueError(f"{column} must be numeric, complete, and finite")
    if plot_data[estimate_col].lt(0).any():
        raise ValueError(f"{estimate_col} must be non-negative")
    if use_intervals:
        if (
            plot_data[ci_low_col].gt(plot_data[estimate_col]).any()
            or plot_data[ci_high_col].lt(plot_data[estimate_col]).any()
        ):
            raise ValueError("Confidence intervals must contain their estimates")

    positions: list[float] = []
    group_centers: list[float] = []
    group_ranges: list[tuple[float, float]] = []
    for group_index in range(len(PAIRED_SBS12_ORDER)):
        start = group_index * (16 + group_gap)
        group_positions = start + np.arange(16, dtype=float)
        positions.extend(group_positions.tolist())
        group_centers.append(float(group_positions.mean()))
        group_ranges.append(
            (float(group_positions[0] - 0.4), float(group_positions[-1] + 0.4))
        )
    x = np.asarray(positions, dtype=float)

    upper_column = ci_high_col if use_intervals else estimate_col
    shared_ymax = float(plot_data[upper_column].max())
    if shared_ymax <= 0:
        shared_ymax = 1.0
    shared_ymax *= 1.15

    n_genes = len(ordered_genes)
    if figsize is None:
        figsize = (24, 3.2 * n_genes + 2.2)
    fig, axes = plt.subplots(
        n_genes,
        1,
        figsize=figsize,
        sharex=True,
        sharey=True,
        squeeze=False,
        constrained_layout=True,
    )
    flat_axes = axes.ravel()
    labels = dict(gene_labels or {})
    bar_colors = [
        SBS12_COLOR_MAPPING[
            parse_context_mutation(context).central_substitution
        ]
        for context in ordered_contexts
    ]

    for axis_index, (axis, gene) in enumerate(zip(flat_axes, ordered_genes)):
        gene_data = (
            plot_data.loc[plot_data[gene_col].eq(gene)]
            .set_index(mutation_col)
            .loc[list(ordered_contexts)]
        )
        estimates = gene_data[estimate_col].to_numpy(dtype=float)
        yerr = None
        if use_intervals:
            yerr = np.vstack(
                (
                    estimates - gene_data[ci_low_col].to_numpy(dtype=float),
                    gene_data[ci_high_col].to_numpy(dtype=float) - estimates,
                )
            )
        axis.bar(
            x,
            estimates,
            width=0.72,
            color=bar_colors,
            alpha=0.9,
            yerr=yerr,
            error_kw={
                "ecolor": "#3F3F46",
                "elinewidth": 0.8,
                "capsize": 1.2,
                "capthick": 0.8,
                "alpha": 0.85,
            },
            zorder=3,
        )
        axis.set_ylim(0, shared_ymax)
        axis.set_xlim(x[0] - 0.8, x[-1] + 0.8)
        axis.grid(False, axis="x")
        axis.grid(axis="y", alpha=0.55, linewidth=0.55, zorder=0)
        axis.set_axisbelow(True)
        axis.set_ylabel(ylabel)

        panel_label = labels.get(gene, str(gene))
        if n_species_col is not None and n_species_col in gene_data.columns:
            sample_sizes = pd.unique(gene_data[n_species_col].dropna())
            if len(sample_sizes) == 1:
                panel_label += f" (N = {int(sample_sizes[0])})"
        axis.set_title(panel_label, fontsize=12, fontweight="bold", loc="left")

        axis_transform = axis.get_xaxis_transform()
        for group_index, substitution in enumerate(PAIRED_SBS12_ORDER):
            start, stop = group_ranges[group_index]
            color = SBS12_COLOR_MAPPING[substitution]
            axis.plot(
                (start, stop),
                (1.015, 1.015),
                color=color,
                linewidth=4.2,
                solid_capstyle="butt",
                transform=axis_transform,
                clip_on=False,
            )
            text_color = "#6B7280" if color == "silver" else color
            axis.text(
                group_centers[group_index],
                1.025,
                substitution,
                color=text_color,
                fontsize=8,
                fontweight="bold",
                ha="center",
                va="bottom",
                transform=axis_transform,
                clip_on=False,
            )
            if group_index < len(PAIRED_SBS12_ORDER) - 1:
                boundary = (stop + group_ranges[group_index + 1][0]) / 2
                axis.axvline(
                    boundary,
                    color="#D1D5DB",
                    linewidth=0.6,
                    zorder=1,
                )

        if axis_index < n_genes - 1:
            axis.tick_params(axis="x", labelbottom=False, bottom=False)

    flat_axes[-1].set_xticks(x)
    flat_axes[-1].set_xticklabels(
        ordered_contexts,
        rotation=90,
        fontsize=ticksize,
        fontname=fontname,
    )
    flat_axes[-1].set_xlabel("Context mutation (heavy-strand orientation)")
    fig.suptitle(title, fontsize=titlesize, fontweight="bold")
    fig.align_ylabels(flat_axes)
    return fig, flat_axes


def summarize_context_tsss_slopes(
    data: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    *,
    genes: Sequence[str],
    central_substitutions: Sequence[str] | None = None,
    contexts: Sequence[str] | None = None,
    tsss_col: str = "dssh_proxy",
    label_col: str | None = "display_gene",
    value_col: str = "MutSpec",
    weighting: str = "species",
    confidence: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | np.random.Generator = 0,
    require_complete_profiles: bool = True,
    minimum_species: int = 2,
    context_order: Sequence[str] = CONTEXT192_ORDER,
    species_col: str = "Species",
    class_col: str = "Class",
    gene_col: str = "Gene",
    mutation_col: str = "Mut",
    expected_col: str = "Expected",
    observed_col: str = "Observed",
) -> pd.DataFrame:
    """Summarize within-species context slopes along a numeric DssH proxy.

    A slope is first calculated within each species across the selected genes.
    Species are then bootstrapped either together or within class.  The
    centered-bootstrap p-value and BH q-value are descriptive screening tools;
    they do not account for phylogeny, linkage among contexts, or uncertainty in
    the fixed DssH proxy.
    """

    ordered_genes = _normalise_genes(genes)
    ordered_contexts = _normalise_context_order(context_order)
    _validate_bootstrap_options("mean", weighting, confidence, n_boot)
    if minimum_species < 2:
        raise ValueError("minimum_species must be at least two")
    if weighting == "equal_class" and class_col not in data.columns:
        raise ValueError(f"weighting='equal_class' requires {class_col!r}")

    if require_complete_profiles:
        validate_matched_context192(
            data,
            genes=ordered_genes,
            context_order=ordered_contexts,
            species_col=species_col,
            class_col=class_col,
            gene_col=gene_col,
            mutation_col=mutation_col,
            value_col=value_col,
            expected_col=expected_col,
            observed_col=observed_col,
        )
    prepared = _prepare_component_rows(
        data,
        genes=ordered_genes,
        value_col=value_col,
        allow_missing_values=not require_complete_profiles,
        species_col=species_col,
        class_col=class_col,
        gene_col=gene_col,
        mutation_col=mutation_col,
    )
    selected_contexts = _select_contexts(
        ordered_contexts,
        central_substitutions=central_substitutions,
        contexts=contexts,
    )
    prepared = prepared.loc[prepared[mutation_col].isin(selected_contexts)]
    if prepared.empty:
        raise ValueError("No rows remain for the selected contexts")

    metadata = _validate_gene_metadata(
        gene_metadata,
        ordered_genes,
        gene_col=gene_col,
        tsss_col=tsss_col,
        label_col=label_col,
    )
    slope_genes = tuple(metadata[gene_col])
    x = metadata[tsss_col].to_numpy(dtype=float)
    centered_x = x - x.mean()
    denominator = float(np.sum(centered_x**2))
    if np.isclose(denominator, 0):
        raise ValueError(f"{tsss_col} must vary across genes")

    rng = _as_rng(random_state)
    alpha = (1 - confidence) / 2
    records: list[dict[str, object]] = []
    for context in selected_contexts:
        component = prepared.loc[prepared[mutation_col].eq(context)]
        if component.empty:
            continue
        matrix = component.pivot(
            index=species_col, columns=gene_col, values=value_col
        ).reindex(columns=list(slope_genes))
        matrix = matrix.dropna(axis=0, how="any")
        if len(matrix) < minimum_species:
            continue

        species_class_rows = component.loc[
            component[species_col].isin(matrix.index),
            [species_col, class_col],
        ].drop_duplicates()
        ambiguity = species_class_rows.groupby(species_col, observed=True)[
            class_col
        ].nunique()
        if ambiguity.gt(1).any():
            examples = list(ambiguity[ambiguity.gt(1)].index[:3])
            raise ValueError(
                f"Species map to multiple {class_col} values; examples: {examples}"
            )
        species_classes = (
            species_class_rows.drop_duplicates(species_col)
            .set_index(species_col)
            .reindex(matrix.index)[class_col]
            .to_numpy()
        )
        slopes = (matrix.to_numpy(dtype=float) @ centered_x) / denominator
        estimate, bootstrap = _bootstrap_group_estimate(
            slopes,
            species_classes,
            estimator=np.mean,
            weighting=weighting,
            n_boot=n_boot,
            rng=rng,
        )
        low, high = np.quantile(bootstrap, [alpha, 1 - alpha])
        centered_bootstrap = bootstrap - estimate
        p_value = (
            1
            + np.count_nonzero(
                np.abs(centered_bootstrap) >= abs(float(estimate))
            )
        ) / (n_boot + 1)
        parsed = parse_context_mutation(context)
        class_sizes = pd.Series(species_classes).value_counts()
        records.append(
            {
                mutation_col: context,
                "Left": parsed.left,
                "Reference": parsed.reference,
                "Alternate": parsed.alternate,
                "Right": parsed.right,
                "CentralSubstitution": parsed.central_substitution,
                "FlankingContext": parsed.flanking_context,
                "ReferenceTrinucleotide": parsed.reference_trinucleotide,
                "n_species": int(len(slopes)),
                "n_classes": int(len(class_sizes)),
                "min_species_per_class": int(class_sizes.min()),
                "max_species_per_class": int(class_sizes.max()),
                "mean_slope_per_proxy_unit": float(estimate),
                "ci_low": float(low),
                "ci_high": float(high),
                "median_species_slope": float(np.median(slopes)),
                "fraction_species_positive": float(np.mean(slopes > 0)),
                "bootstrap_p_two_sided": float(p_value),
                "confidence": confidence,
                "weighting": weighting,
                "value_col": value_col,
                "p_value_method": "centered_species_bootstrap",
            }
        )

    result = pd.DataFrame.from_records(records)
    if result.empty:
        raise ValueError(
            "No selected context has enough species with values for every gene"
        )
    result["fdr_bh_q"] = _benjamini_hochberg(
        result["bootstrap_p_two_sided"].to_numpy(dtype=float)
    )
    return result


def _normalise_genes(genes: Sequence[str]) -> tuple[str, ...]:
    ordered = tuple(genes)
    if len(ordered) < 2:
        raise ValueError("genes must contain at least two entries")
    if len(set(ordered)) != len(ordered):
        raise ValueError("genes must contain unique entries")
    return ordered


def _normalise_context_order(context_order: Sequence[str]) -> tuple[str, ...]:
    contexts = tuple(str(item) for item in context_order)
    if len(contexts) != 192 or len(set(contexts)) != 192:
        raise ValueError("context_order must contain exactly 192 unique labels")
    for context in contexts:
        parse_context_mutation(context)
    if set(contexts) != set(CONTEXT192_ORDER):
        missing = sorted(set(CONTEXT192_ORDER) - set(contexts))[:5]
        extra = sorted(set(contexts) - set(CONTEXT192_ORDER))[:5]
        raise ValueError(
            "context_order must contain the canonical 192 components; "
            f"missing={missing}, extra={extra}"
        )
    return contexts


def _select_contexts(
    context_order: Sequence[str],
    *,
    central_substitutions: Sequence[str] | None,
    contexts: Sequence[str] | None,
) -> tuple[str, ...]:
    selected = list(context_order)
    if central_substitutions is not None:
        substitutions = tuple(str(item) for item in central_substitutions)
        if not substitutions or len(set(substitutions)) != len(substitutions):
            raise ValueError(
                "central_substitutions must contain at least one unique label"
            )
        unknown = set(substitutions).difference(SBS12_ORDER)
        if unknown:
            raise ValueError(f"Unknown central substitutions: {sorted(unknown)}")
        substitution_set = set(substitutions)
        selected = [
            context
            for context in selected
            if parse_context_mutation(context).central_substitution
            in substitution_set
        ]
    if contexts is not None:
        requested = tuple(str(item) for item in contexts)
        if not requested or len(set(requested)) != len(requested):
            raise ValueError("contexts must contain at least one unique label")
        unknown = set(requested).difference(context_order)
        if unknown:
            raise ValueError(f"Unknown contexts: {sorted(unknown)[:5]}")
        requested_set = set(requested)
        selected = [item for item in selected if item in requested_set]
    if not selected:
        raise ValueError("No contexts satisfy the requested filters")
    return tuple(selected)


def _prepare_component_rows(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    value_col: str,
    allow_missing_values: bool,
    species_col: str,
    class_col: str,
    gene_col: str,
    mutation_col: str,
) -> pd.DataFrame:
    required = {species_col, class_col, gene_col, mutation_col, value_col}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    subset = data.loc[data[gene_col].isin(genes)].copy()
    absent = [gene for gene in genes if gene not in set(subset[gene_col])]
    if absent:
        raise ValueError(f"Genes absent from the analysis rows: {absent}")
    if subset[[species_col, class_col, gene_col, mutation_col]].isna().any().any():
        raise ValueError("Required component identifiers must not be missing")
    duplicate = subset.duplicated([species_col, gene_col, mutation_col])
    if duplicate.any():
        examples = subset.loc[
            duplicate, [species_col, gene_col, mutation_col]
        ].head(3)
        raise ValueError(
            "Duplicate Species/Gene/Mut rows were found; examples: "
            f"{examples.to_dict('records')}"
        )

    numeric = pd.to_numeric(subset[value_col], errors="coerce")
    invalid_nonmissing = subset[value_col].notna() & numeric.isna()
    if invalid_nonmissing.any():
        raise ValueError(f"Non-missing {value_col} values must be numeric")
    if not np.isfinite(numeric.dropna()).all():
        raise ValueError(f"Non-missing {value_col} values must be finite")
    if numeric.dropna().lt(0).any():
        raise ValueError(f"Non-missing {value_col} values must be non-negative")
    if numeric.isna().any():
        if not allow_missing_values:
            raise ValueError(f"{value_col} must be complete")
        subset = subset.loc[numeric.notna()].copy()
        numeric = numeric.loc[numeric.notna()]
    subset[value_col] = numeric.to_numpy(dtype=float)

    for label in pd.unique(subset[mutation_col]):
        parse_context_mutation(str(label))
    species_classes = subset[[species_col, class_col]].drop_duplicates()
    ambiguity = species_classes.groupby(species_col, observed=True)[
        class_col
    ].nunique()
    if ambiguity.gt(1).any():
        examples = list(ambiguity[ambiguity.gt(1)].index[:3])
        raise ValueError(
            f"Species map to multiple {class_col} values; examples: {examples}"
        )
    return subset


def _validate_bootstrap_options(
    estimator: str,
    weighting: str,
    confidence: float,
    n_boot: int,
) -> None:
    if estimator not in {"mean", "median"}:
        raise ValueError("estimator must be 'mean' or 'median'")
    if weighting not in {"species", "equal_class"}:
        raise ValueError("weighting must be 'species' or 'equal_class'")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    if n_boot < 1:
        raise ValueError("n_boot must be at least one")


def _bootstrap_matched_profile_estimate(
    data: pd.DataFrame,
    *,
    genes: Sequence[str],
    contexts: Sequence[str],
    value_col: str,
    estimator,
    weighting: str,
    n_boot: int,
    rng: np.random.Generator,
    species_col: str,
    class_col: str,
    gene_col: str,
    mutation_col: str,
    chunk_size: int = 50,
) -> tuple[tuple[object, ...], np.ndarray, np.ndarray, np.ndarray]:
    """Bootstrap matched species as joint gene-by-context profiles.

    The same sampled species indices are applied to every gene and context in
    one replicate.  The returned bootstrap array has shape
    ``(n_boot, n_genes, n_contexts)``.
    """

    ordered_genes = tuple(genes)
    ordered_contexts = tuple(contexts)
    species = tuple(sorted(pd.unique(data[species_col]), key=str))
    if not species:
        raise ValueError("Cannot bootstrap an empty matched species cohort")

    species_class_rows = data[[species_col, class_col]].drop_duplicates()
    ambiguity = species_class_rows.groupby(species_col, observed=True)[
        class_col
    ].nunique()
    if ambiguity.gt(1).any():
        examples = list(ambiguity[ambiguity.gt(1)].index[:3])
        raise ValueError(
            f"Species map to multiple {class_col} values; examples: {examples}"
        )
    species_classes = (
        species_class_rows.drop_duplicates(species_col)
        .set_index(species_col)
        .reindex(species)[class_col]
    )
    if species_classes.isna().any():
        raise ValueError("Every matched species must have a class label")
    class_values = species_classes.to_numpy()

    cube = np.empty(
        (len(species), len(ordered_genes), len(ordered_contexts)),
        dtype=float,
    )
    for gene_index, gene in enumerate(ordered_genes):
        matrix = (
            data.loc[data[gene_col].eq(gene)]
            .pivot(index=species_col, columns=mutation_col, values=value_col)
            .reindex(index=species, columns=ordered_contexts)
        )
        if matrix.isna().any().any():
            raise ValueError(
                "Joint profile bootstrap requires every matched "
                f"species/gene/context row; incomplete gene={gene!r}"
            )
        cube[:, gene_index, :] = matrix.to_numpy(dtype=float)

    bootstrap = np.empty(
        (n_boot, len(ordered_genes), len(ordered_contexts)), dtype=float
    )
    if weighting == "species":
        estimate = estimator(cube, axis=0)
        for start in range(0, n_boot, chunk_size):
            stop = min(start + chunk_size, n_boot)
            draws = rng.integers(
                0, len(species), size=(stop - start, len(species))
            )
            bootstrap[start:stop] = estimator(cube[draws], axis=1)
    else:
        class_indices = [
            np.flatnonzero(class_values == item)
            for item in sorted(pd.unique(class_values), key=str)
        ]
        estimate = np.mean(
            [estimator(cube[index], axis=0) for index in class_indices],
            axis=0,
        )
        for start in range(0, n_boot, chunk_size):
            stop = min(start + chunk_size, n_boot)
            class_bootstrap = []
            for index in class_indices:
                draws = rng.integers(
                    0, len(index), size=(stop - start, len(index))
                )
                class_bootstrap.append(estimator(cube[index][draws], axis=1))
            bootstrap[start:stop] = np.mean(class_bootstrap, axis=0)

    if (
        estimator is np.mean
        and len(ordered_contexts) == 192
        and set(ordered_contexts) == set(CONTEXT192_ORDER)
    ):
        if not np.allclose(estimate.sum(axis=1), 1.0, rtol=0, atol=1e-8):
            raise ValueError("Mean matched profile estimates must sum to one")
        if not np.allclose(
            bootstrap.sum(axis=2), 1.0, rtol=0, atol=1e-8
        ):
            raise ValueError(
                "Every mean matched-profile bootstrap replicate must sum to one"
            )

    return species, class_values, np.asarray(estimate), bootstrap


def _bootstrap_group_estimate(
    values: np.ndarray,
    classes: np.ndarray,
    *,
    estimator,
    weighting: str,
    n_boot: int,
    rng: np.random.Generator,
    chunk_size: int = 250,
) -> tuple[float, np.ndarray]:
    values = np.asarray(values, dtype=float)
    classes = np.asarray(classes)
    if len(values) == 0:
        raise ValueError("Cannot bootstrap an empty group")
    bootstrap = np.empty(n_boot, dtype=float)

    if weighting == "species":
        estimate = float(estimator(values))
        for start in range(0, n_boot, chunk_size):
            stop = min(start + chunk_size, n_boot)
            draws = rng.integers(
                0, len(values), size=(stop - start, len(values))
            )
            bootstrap[start:stop] = estimator(values[draws], axis=1)
        return estimate, bootstrap

    class_values = [
        values[classes == item] for item in sorted(pd.unique(classes), key=str)
    ]
    if not class_values:
        raise ValueError("Equal-class weighting requires at least one class")
    estimate = float(np.mean([estimator(group) for group in class_values]))
    for start in range(0, n_boot, chunk_size):
        stop = min(start + chunk_size, n_boot)
        class_bootstrap = np.empty(
            (stop - start, len(class_values)), dtype=float
        )
        for class_index, group in enumerate(class_values):
            draws = rng.integers(
                0, len(group), size=(stop - start, len(group))
            )
            class_bootstrap[:, class_index] = estimator(
                group[draws], axis=1
            )
        bootstrap[start:stop] = class_bootstrap.mean(axis=1)
    return estimate, bootstrap


def _validate_gene_metadata(
    gene_metadata: pd.DataFrame,
    genes: Sequence[str],
    *,
    gene_col: str,
    tsss_col: str,
    label_col: str | None,
) -> pd.DataFrame:
    required = {gene_col, tsss_col}
    if label_col is not None:
        required.add(label_col)
    missing = required.difference(gene_metadata.columns)
    if missing:
        raise ValueError(f"gene_metadata is missing columns: {sorted(missing)}")
    columns = [gene_col, tsss_col]
    if label_col is not None:
        columns.append(label_col)
    metadata = gene_metadata.loc[
        gene_metadata[gene_col].isin(genes), columns
    ].copy()
    if metadata[gene_col].duplicated().any():
        raise ValueError("gene_metadata must contain one row per gene")
    if set(metadata[gene_col]) != set(genes):
        absent = sorted(set(genes) - set(metadata[gene_col]), key=str)
        raise ValueError(f"gene_metadata is missing genes: {absent}")
    metadata[tsss_col] = pd.to_numeric(metadata[tsss_col], errors="coerce")
    if metadata[tsss_col].isna().any() or not np.isfinite(metadata[tsss_col]).all():
        raise ValueError(f"{tsss_col} values must be numeric, complete, and finite")
    if metadata[tsss_col].duplicated().any():
        raise ValueError(f"{tsss_col} values must be unique across genes")
    if label_col is not None and metadata[label_col].isna().any():
        raise ValueError(f"{label_col} values must be complete")
    return metadata.sort_values(tsss_col).reset_index(drop=True)


def _as_rng(
    random_state: int | np.random.Generator,
) -> np.random.Generator:
    if isinstance(random_state, np.random.Generator):
        return random_state
    return np.random.default_rng(random_state)


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    if p_values.ndim != 1 or not np.isfinite(p_values).all():
        raise ValueError("p-values must be a finite one-dimensional array")
    if ((p_values < 0) | (p_values > 1)).any():
        raise ValueError("p-values must lie between zero and one")
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return result


__all__ = [
    "BASE_ORDER",
    "CONTEXT192_ORDER",
    "PAIRED_CONTEXT192_ORDER",
    "PAIRED_SBS12_ORDER",
    "SBS12_COLOR_MAPPING",
    "Context192ProfileInfo",
    "ContextMutation",
    "QC_STATUS_ORDER",
    "SBS12_ORDER",
    "add_context_columns",
    "annotate_context_qc",
    "calculate_within_substitution_shares",
    "orient_contexts_to_heavy_strand",
    "parse_context_mutation",
    "plot_context192_bar_spectra",
    "reverse_complement_context",
    "select_common_species_192",
    "summarize_context_components",
    "summarize_context_tsss_slopes",
    "validate_matched_context192",
]
