# TBSS-Signature

This project studies whether mitochondrial mutation spectra vary with the time
that mtDNA remains single-stranded. The repository is organised as a sequence
of self-contained analyses. Each numbered folder contains its own notebook or
script, generated tables, figures, and a README that documents that stage.

## Project structure

```text
1init_data/
|-- data/                         source datasets
|-- Check192spec.ipynb            input quality-control notebook
`-- README.md
2species_intersection/
|-- data/                         intersection tables
|-- figures/                      intersection plots
|-- count_common_species.py       reusable intersection code
|-- CommonSpeciesStatistics.ipynb
`-- README.md
3compare_t_genes/
|-- data/                         matched cohorts and summaries
|-- figures/                      spectrum comparison plots
|-- mutation_comparison.py        reusable analysis and plotting functions
|-- CompareTGenes.ipynb
`-- README.md
4tsss_gradient/
|-- data/                         TSSS metadata and summaries
|-- figures/                      mutation-gradient plots
|-- TSSSMutationGradient.ipynb
`-- README.md
```

## Analysis order

1. `1init_data` documents and checks the supplied mutation-spectrum tables.
2. `2species_intersection` counts the mammals shared by every gene combination
   and records the exact members of each matched cohort.
3. `3compare_t_genes` compares 12-component spectra for two to five genes using
   the same species for every selected gene.
4. `4tsss_gradient` explores changes in selected mutation types along an
   explicit proxy for single-stranded duration.

The main input for stages 2--4 is
`1init_data/data/MutSpecVertebrates12.csv.gz`. Run notebooks in the numbered
order when regenerating the complete analysis. Tables are written to the
corresponding `data/` folder and plots to `figures/`; source data remain in
`1init_data/data/`.

See the README inside each stage for its assumptions, commands, outputs, and
interpretation limits.
