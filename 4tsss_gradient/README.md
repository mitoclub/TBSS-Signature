# 4. Mutation gradient along single-stranded duration

## Purpose

This analysis visualises how selected mutation-spectrum components vary among
mitochondrial genes positioned at different distances from the light-strand
origin. It is an exploratory within-species description of a proposed
single-stranded-duration gradient.

## Inputs

- Spectra: `../1init_data/data/MutSpecVertebrates12.csv.gz`.
- Reusable validation, orientation, plotting, and bootstrap functions:
  `../3compare_t_genes/mutation_comparison.py`.

The current analysis uses the exact shared mammal cohort for `CO1`, `CO3`, and
`Cytb`. It examines `C>T`, `A>G`, `G>A`, and `T>C` after orienting substitutions
to heavy-strand notation.

## TSSS proxy

`TSSSMutationGradient.ipynb` records rCRS start and end positions for the three
genes, calculates each midpoint, and uses

```text
DssH proxy = 2 * (gene midpoint - OriL start) / mtDNA length
```

with an mtDNA length of 16,569 bp and an OriL start of 5,730. The explicit
numeric metadata are saved with the results; gene order alone is never used as
the explanatory variable.

## Workflow

The notebook:

1. constructs and saves the gene-level TSSS metadata;
2. selects species present for all three genes;
3. validates the matched 12-component profiles;
4. plots species trajectories plus bootstrap means and intervals; and
5. summarizes descriptive within-species slopes along the proxy.

Run `4tsss_gradient/TSSSMutationGradient.ipynb` after the earlier stages. To
change the genes or mutation types, edit `GENES`, `MUTATIONS_TO_PLOT`, and the
corresponding metadata table together.

## Outputs

The `data/` folder contains:

- `gene_tsss_proxy.csv`: positions, labels, and calculated proxy values;
- `matched_species_CO1_CO3_Cytb.csv`: exact three-gene cohort;
- `tsss_gradient_gene_summary.csv`: gene/mutation estimates and intervals;
- `tsss_gradient_slope_summary.csv`: descriptive slope summaries.

The `figures/` folder contains `tsss_gradient_selected_mutations.png`.

## Caveats

- The proxy is derived from fixed human rCRS coordinates and is not a
  species-specific replication-timing measurement.
- With only three genes, slope estimates cannot separate TSSS from gene
  identity, base composition, selection, or other genomic factors.
- The bootstrap preserves matched species but does not model phylogenetic
  relatedness. Results should therefore be interpreted as descriptive, not
  causal.
