# 3. Compare mutation spectra between genes

## Purpose

This analysis compares context-free 12-component mutation spectra across two to
five mitochondrial genes while holding the species cohort constant. The current
project scope is `Mammalia`. The notebook example uses `CO1` and `Cytb`, the
available requested pair with the largest matched mammalian cohort: 104 species.

## Inputs

- Primary spectra: `../1init_data/data/MutSpecVertebrates12.csv.gz`.
- Species-overlap tables from `../2species_intersection/data/` document which
  gene combinations have adequate matched cohorts.

`CompareTGenes.ipynb` currently selects the exact common species directly from
the primary table and saves that cohort alongside the results.

## Workflow

The notebook selects `Mammalia`, forms the all-gene intersection, orients
substitution labels to heavy-strand notation, validates
every matched profile, and then draws:

- the complete 12-component spectrum comparison;
- detailed paired comparisons for `C>T` and `G>A`; and
- paired boxplots for the within-profile ratios `C>T/G>A` and `A>G/T>C`.

The ratios match stage 4: they are calculated separately within every
species/gene profile, without a pseudocount, rather than as ratios of
cross-species means. Because profile normalization cancels in a ratio, the
values are identical for normalized and pre-normalization spectrum components.

The boxplots use a clean, colorblind-safe journal style. Every faint trajectory
is retained and represents one matched species. Ratio panels use a logarithmic
y-axis, with a reference line at ratio 1, so their long right tails do not
compress the boxes. Exact zero numerators are kept at zero in a short linear
segment attached to the logarithmic scale; no plotting pseudocount is added.

The upper-right annotation reports a paired Wilcoxon signed-rank test because
the same species are measured for every gene; Mann--Whitney U would be an
unpaired test and is therefore not used. Directions are specified before the
test and describe Gene2 relative to Gene1: the current notebook tests
`CytB > COX1` for both ratios and `C>T`, and `CytB < COX1` for `G>A`. With more
than two selected genes, all gene pairs are tested and displayed p-values use
Holm correction within each plotted quantity. A numerical underflow is shown as
an inequality rather than `p = 0`.

Means and uncertainty intervals are calculated by resampling whole species, so
the matched relationship between genes and mutation components is preserved.

To analyse another combination, change `GENES` in `CompareTGenes.ipynb` to any
two to five available genes. The notebook will regenerate filenames from that
selection. `CLASS_FILTER` remains fixed to `"Mammalia"` for the current project.

## Reusable functions

`mutation_comparison.py` provides strict reusable helpers including
`select_common_species`, `orient_substitutions_to_heavy_strand`,
`validate_matched_spectra`, `plot_matched_spectra`, and
`plot_mutation_comparison`. It also provides
`calculate_mutation_ratios`, `plot_mutation_ratio_comparison`,
`plot_tsss_ratio_gradient`, and
`summarize_tsss_ratio_slopes` for ratios calculated separately within each
species/gene profile. The default ratios are `C>T/G>A` and `A>G/T>C`; a zero
denominator raises an error and no pseudocount is added silently.

Both paired boxplot functions accept `pairwise_alternative`. Its compatibility
default is `"two-sided"`; `"greater"` means Gene2 > Gene1 and `"less"` means
Gene2 < Gene1. The returned test table records this public direction, SciPy's
internally reversed alternative on Gene1 - Gene2, and the definition of the
reported paired difference.

The summary and TSSS-gradient helpers accept `require_sum_to_one=False` for a
non-compositional value column. Stage 4 uses this explicitly for
`Observed.fillna(0) / Expected`; callers must supply an honest axis label because
that quantity is an opportunity-adjusted reconstructed burden, not an absolute
mutation rate.

The plotting functions expect a DataFrame already restricted to exactly the
same species for every selected gene. They validate that cohort but do not
silently intersect or drop species. From the repository root, add
`3compare_t_genes/` to `sys.path` before importing `mutation_comparison`.

## Outputs

Current reproducible outputs in `data/` are:

- `common_species_<genes>.csv`: identifiers and source classes in the matched
  cohort;
- `matched_spectra_<genes>.csv`: row-level matched spectra;
- `matched_12_component_summary_<genes>.csv`: estimates and bootstrap
  confidence intervals;
- `matched_mutation_ratios_<genes>.csv`: per-species, per-gene ratio values;
- `matched_mutation_ratio_summary_<genes>.csv`: ratio means, medians, and
  species-bootstrap confidence intervals;
- `matched_CT_pairwise_tests_<genes>.csv`: paired test for the detailed `C>T`
  boxplot;
- `matched_GA_pairwise_tests_<genes>.csv`: paired test for the detailed `G>A`
  boxplot;
- `matched_mutation_ratio_pairwise_tests_<genes>.csv`: paired tests for the
  ratio boxplots;
- `candidate_pair_counts.csv`: availability of the requested candidate pairs
  in the current taxonomic scope.

Current figures in `figures/` are `matched_12_component_<genes>.png` plus PNG
and vector PDF versions of `matched_CT_<genes>`, `matched_GA_<genes>`, and
`matched_mutation_ratios_<genes>`. Earlier mammal-only files with `paired_` or
`mean_` prefixes are preserved under `legacy_mammalia/` for provenance and are
not regenerated by the current notebook.

## Caveats

- `ND6` is absent from the supplied 12-component table and is not substituted
  with `ND2`.
- Heavy-strand complementation must be applied exactly once; applying it to an
  already oriented table reverses the intended notation.
- Bootstrap intervals describe variation across the matched species. They do
  not correct for phylogenetic non-independence.
- Ratios can be skewed and sensitive to small denominators; the saved
  species-level values should be inspected together with means and medians.
- The Wilcoxon tests use species as paired observations but do not account for
  phylogenetic non-independence; p-values should therefore be interpreted as
  descriptive evidence for this supplied cohort.
- The 104-species cohort is mammalian but remains taxonomically and
  phylogenetically non-independent.
