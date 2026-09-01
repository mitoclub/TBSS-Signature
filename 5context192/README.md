# 5. Matched-gene analysis of 192-component spectra

## Purpose

This stage asks whether the trinucleotide context of four heavy-strand
transitions (`C>T`, `A>G`, `G>A`, and `T>C`) changes across the matched
`CO1`--`CO3`--`Cytb` gene series. It also provides a compact quality-control
view of the complete 192-component spectra.

The current scope is `Mammalia`. The exact intersection contains 52 species
with a complete spectrum for every one of the three genes; comparisons are
therefore within the same species cohort.

## Inputs

- `../1init_data/data/MutSpecVertebrates192.csv.gz`: probabilistically weighted
  reconstructed substitutions, expected synonymous opportunities, and
  normalized 192-component spectra.
- `context192_analysis.py`: reusable parsing, heavy-strand orientation,
  validation, context-share, bootstrap-summary, and positional-slope helpers.

Run `Context192MatchedGenes.ipynb` from either the repository root or this
folder. The notebook creates `data/` and `figures/` automatically.

## Orientation and metrics

Source context labels have the form `L[reference>alternate]R`. The notebook
applies one full reverse-complement transformation to the complete label:
flanks are swapped and complemented, and both central alleles are
complemented. The original label is retained as `Mut_raw`, and the oriented
table is explicitly marked so that a second orientation call is rejected.

The complete-spectrum figure follows the reference repository's
`plot_mutspec192` layout. Its 12 blocks are ordered
`C>A, G>T, C>G, G>C, C>T, G>A, T>A, A>T, T>C, A>G, T>G, A>C`; each block
contains 16 contexts, complementary blocks use reverse-complement-aligned
context order, and each complementary pair shares a color. The legacy
PyMutSpec package is not imported because its dependency chain is incompatible
with the current Python version; `context192_analysis.py` reproduces the
ordering and plotting logic directly.

Three related quantities are kept distinct:

1. **Whole-spectrum weight (`MutSpec`)** is the supplied normalized component;
   all 192 components in a species/gene profile sum to one.
2. **Within-substitution context share** divides a component by the sum of the
   16 contexts having the same central substitution in the same species and
   gene. These 16 shares sum to one when that substitution has non-zero total
   weight. This view describes context allocation, not absolute mutation
   burden. A zero-total substitution remains undefined (no pseudocount is
   added), so the conditional summaries retain component-specific sample
   sizes.
3. **Observed/expected sensitivity** is `Observed / Expected` only where
   `Expected > 0`, with missing reconstructed weights treated as zero. It is an
   opportunity-normalized descriptive burden, not a per-generation
   mutation rate. It is used only as a sensitivity table and is never imputed
   for rows with no expected opportunity.

The 192-component table and the 12-component table were opportunity-corrected
and normalized independently. Consequently, summing 16 normalized context
weights is a marginal of the 192-component model, not an exact reconstruction
of the supplied 12-component `MutSpec`.

## Quality control

Each valid row is assigned to one of three mutually exclusive states:

- `no_expected_opportunity`: no expected synonymous opportunity for that
  context;
- `opportunity_observed_zero`: an opportunity exists but no mutation was
  observed;
- `opportunity_observed_positive`: both opportunity and observation are
  positive.

The helper also reserves
`no_expected_opportunity_but_observed_positive` as an invalid fourth state;
strict validation stops if it occurs. A `no_expected_opportunity` row is an
opportunity/identifiability zero in this coding
sequence analysis; it must not be interpreted as evidence that the molecular
mutation cannot occur. The notebook checks that every matched gene/species
profile contains exactly 192 unique labels, has finite non-negative `MutSpec`
values, and sums to one.

## Workflow

The notebook:

1. filters to `Mammalia` and saves the exact 52-species, three-gene cohort;
2. validates the raw profiles, reverse-complements contexts once, and records
   the opportunity/observation QC states;
3. plots the mean complete spectra for all three genes as three aligned
   reference-style 192-bar panels with species-bootstrap intervals;
4. plots absolute whole-spectrum weights and within-substitution shares for the
   four selected transitions across all 16 flanking contexts;
5. estimates matched-species context slopes for both whole-spectrum weights
   and within-substitution shares along the same fixed human-rCRS DssH proxy
   used in stage 4; and
6. reports an Expected-positive observed/expected sensitivity analysis.

For complete spectra, bootstrap intervals resample matched species as whole
profiles: one draw is reused across all three genes and all 192 components in
each replicate. Consequently, mean-spectrum bootstrap replicates remain
normalized to one and preserve the matched covariance needed for later gene
contrasts. In the deliberately incomplete `Observed / Expected` sensitivity,
available species differ by component, so those intervals are calculated
component by component and report component-specific sample sizes. The slope
field `fdr_bh_q`
applies Benjamini--Hochberg false-discovery-rate correction to descriptive
centered-bootstrap p-values across the 64 selected context tests. These remain
exploratory screening summaries rather than confirmatory phylogenetic
inference.

## Outputs

The `data/` folder contains:

- `gene_tsss_proxy.csv`: fixed rCRS coordinates and calculated proxy values;
- `matched_species_CO1_CO3_Cytb.csv` and
  `matched_class_counts_CO1_CO3_Cytb.csv`: cohort membership and composition;
- `matched_192_spectra_CO1_CO3_Cytb.csv.gz`: the validated, once-oriented
  matched profiles with QC annotations;
- `context192_qc_summary.csv`: counts of the three QC states;
- `context192_mean_spectrum.csv`: mean complete normalized spectra;
- `transition_context_mutspec_summary.csv`: absolute normalized component
  summaries for the four transitions;
- `transition_context_share_summary.csv`: within-substitution context shares;
- `transition_context_oe_sensitivity.csv`: `Observed / Expected` summaries
  restricted to `Expected > 0`;
- `context_tsss_slopes_species_weighted.csv`: mammalian context gradients;
- `context_tsss_slopes_within_substitution_share.csv`: gradients after
  conditioning on the total weight of each central substitution; and
- `context_tsss_slope_sensitivity.csv`: the three views side by side.

The `figures/` folder contains:

- `mean_192_component_spectra.png`: reference-style 192-bar spectra for COX1,
  COX3, and CytB on a common y-scale;
- `transition_context_mutspec_heatmaps.png`;
- `transition_context_share_heatmaps.png`; and
- `transition_context_tsss_slopes.png`.

## COSMIC signature assignment substage

The [`signatures/`](signatures/) subfolder applies the reference
`4signatures` workflow to each of the three mean matched-gene spectra with two
independent implementations: SigProfilerAssignment in Python and mSigAct in R.
The strand-specific SBS192 profiles are converted to canonical SBS96 only at
the COSMIC fitting boundary. It prepares the same high-transition,
low-transition, and positive high-minus-low variants, with and without
transversions, and retains all 18 gene/variant assignments. See the subfolder
README for exact priors, exclusions, versions, outputs, and interpretation
limits.

## Interpretation limits

- The x-axis uses fixed human rCRS gene coordinates and OriL for every sampled
  vertebrate. It is a positional proxy, not species-specific time spent
  single-stranded.
- The 52 mammalian species remain phylogenetically non-independent and unevenly
  distributed among mammalian clades.
- Only three genes define each slope, so gene identity, sequence composition,
  selection, position, and replication exposure remain confounded.
- `MutSpec` is compositional. A higher component weight can reflect change in
  that component, change elsewhere in the spectrum, or both.
- Some contexts have no synonymous expected opportunity under the genetic code
  and observed codons. Such rows are excluded, not assigned a numerical value,
  in the observed/expected sensitivity analysis.
- Many other contexts have a small positive `Expected` value, so their
  `Observed / Expected` burden can be unstable. That sparse analysis is kept as
  a sensitivity table with component-specific sample sizes, not used as the
  primary context-gradient result.
