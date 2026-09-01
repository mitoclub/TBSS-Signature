# 4. Mutation gradient along a positional DssH proxy

## Purpose

This analysis visualises how selected mutation-spectrum components vary among
mitochondrial genes positioned at different distances from the light-strand
origin. It compares three related views of the same matched species: normalized
12-component spectrum weights, pre-normalization opportunity-adjusted burdens,
and complementary-substitution ratios. It is an exploratory within-species
description of a proposed single-stranded-duration gradient.

## Inputs

- Spectra: `../1init_data/data/MutSpecVertebrates12.csv.gz`.
- Reusable validation, orientation, plotting, and bootstrap functions:
  `../3compare_t_genes/mutation_comparison.py`.

The current analysis uses the exact mammalian shared cohort for `CO1`, `CO3`,
and `Cytb`. It contains 52 species and examines `C>T`, `A>G`, `G>A`, and `T>C`
after orienting substitutions to heavy-strand notation. `CLASS_FILTER` is fixed
to `"Mammalia"` for the current project.

## TSSS proxy

`TSSSMutationGradient.ipynb` records rCRS start and end positions for the three
genes, calculates each midpoint, and uses

```text
DssH proxy = 2 * (gene midpoint - OriL start) / mtDNA length
```

with an mtDNA length of 16,569 bp and an OriL start of 5,730. The explicit
numeric metadata are saved with the results; gene order alone is never used as
the explanatory variable. These are fixed human rCRS coordinates applied as a
common positional scale, not observed single-stranded durations for the sampled
species.

## Spectrum quantities

The normalized plots use the supplied `MutSpec` weights, which sum to one within
each species-gene profile. The unnormalized variation uses the intermediate
quantity from the source derivation:

```text
OpportunityAdjusted = Observed.fillna(0) / Expected
```

This removes the final 100% normalization and adjusts reconstructed substitutions
for synonymous opportunities. It is best described as an
**opportunity-adjusted reconstructed substitution burden**. It is not an
absolute mutation frequency, probability, or per-generation rate. The notebook
explicitly requires finite, strictly positive `Expected` values.

The two ratios are calculated separately within every species and gene:

```text
C>T / G>A
A>G / T>C
```

No pseudocount is used; the analysis stops if a denominator is non-positive.
Because the same normalization constant multiplies every component of a profile,
the ratio is identical for normalized and pre-normalization values:
`(r_i / Z) / (r_j / Z) = r_i / r_j`. The notebook checks this identity
numerically rather than taking a ratio of cross-species means.

## Workflow

The notebook:

1. constructs and saves the gene-level TSSS metadata;
2. selects species present for all three genes;
3. validates the matched 12-component profiles;
4. plots normalized species trajectories plus bootstrap means and intervals;
5. repeats the four-substitution gradient using `Observed / Expected` without
   100% normalization;
6. plots the two per-species complementary-substitution ratios; and
7. summarizes descriptive within-species slopes for all three views.

Run `4tsss_gradient/TSSSMutationGradient.ipynb` after the earlier stages. To
change the genes or mutation types, edit `GENES`, `MUTATIONS_TO_PLOT`, and the
corresponding metadata table together.

## Outputs

The `data/` folder contains:

- `gene_tsss_proxy.csv`: positions, labels, and calculated proxy values;
- `matched_species_CO1_CO3_Cytb.csv`: exact three-gene cohort;
- `matched_class_counts_CO1_CO3_Cytb.csv`: class composition of that cohort;
- `tsss_gradient_gene_summary.csv`: gene/mutation estimates and intervals;
- `tsss_gradient_slope_summary.csv`: normalized-spectrum slope summaries;
- `tsss_opportunity_adjusted_species_values.csv`: row-level `Observed`,
  `Expected`, and `Observed / Expected` values for the four substitutions;
- `tsss_opportunity_adjusted_gene_summary.csv`: unnormalized gene/mutation
  estimates and intervals;
- `tsss_opportunity_adjusted_slope_summary.csv`: unnormalized slope summaries;
- `tsss_mutation_ratio_species_values.csv`: both ratios and their numerator and
  denominator values for each matched species and gene;
- `tsss_mutation_ratio_gene_summary.csv`: ratio estimates and intervals; and
- `tsss_mutation_ratio_slope_summary.csv`: ratio slope summaries.

The `figures/` folder contains:

- `tsss_gradient_selected_mutations.png`: normalized spectrum weights;
- `tsss_opportunity_adjusted_selected_mutations.png`: four
  opportunity-adjusted burdens without 100% normalization; and
- `tsss_mutation_ratio_gradient.png`: the two within-profile ratios.

## Caveats

- The proxy is derived from fixed human rCRS coordinates and is not a
  species-specific replication-timing measurement. Applying the same human
  coordinates to all sampled mammals is still a simplifying assumption;
  differences in mitogenome organisation are not represented.
- With only three genes, slope estimates cannot separate TSSS from gene
  identity, base composition, selection, or other genomic factors.
- Removing the 100% normalization avoids compositional closure but does not turn
  the reconstructed `Observed / Expected` quantity into an absolute mutation
  frequency or rate.
- Ratios can be skewed and sensitive to small denominators. Species-level
  trajectories and denominator values are retained, and arithmetic means should
  be read together with the medians and slope summaries.
- The bootstrap preserves matched species but does not model phylogenetic
  relatedness. Results should therefore be interpreted as descriptive, not
  causal.
