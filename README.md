# TBSS-Signature

This project studies whether mitochondrial mutation spectra vary along a
positional proxy for the time that mtDNA remains single-stranded. The current
taxonomic scope is `Mammalia` throughout the project. The supplied mutation
tables retain other vertebrate classes as immutable source data, but every
analysis filters to mammals before forming cohorts. The repository is organised
as a sequence of self-contained analyses. Each numbered folder
contains its own notebook or script, generated tables, figures, and a README
that documents that stage.

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
5context192/
|-- data/                         matched context spectra and summaries
|-- figures/                      192-component and context-gradient plots
|-- signatures/                   Python/R COSMIC signature assignment
|-- context192_analysis.py        reusable context-spectrum functions
|-- Context192MatchedGenes.ipynb
`-- README.md
6amino_acid_shift/
|-- data/midori/                  cached MIDORI2 GB272 FASTA archives
|-- data/derived_midori/          source manifest and final QC sequences
|-- figures/                      ratio distributions and paired effects
|-- midori_analysis.py            Mammalia discovery, pairing, and selection
|-- amino_acid_shift.py           translation-table and sequence-QC helpers
|-- run_analysis.py               MIDORI extraction/QC entry point
|-- AminoAcidShift.ipynb
`-- README.md
```

## Analysis order

1. `1init_data` documents and checks the supplied mutation-spectrum tables.
2. `2species_intersection` counts mammalian species shared by every gene
   combination and records the exact members of each matched cohort.
3. `3compare_t_genes` compares 12-component spectra for two to five genes using
   the same species for every selected gene. The current `CO1`/`Cytb` example
   contains 104 mammalian species and includes paired
   boxplots for `C>T`, `G>A`, and the `C>T/G>A` and `A>G/T>C` within-profile
   ratios, with pre-specified one-sided paired Wilcoxon tests.
4. `4tsss_gradient` explores changes in selected mutation types along an
   explicit proxy for single-stranded duration. Its current three-gene cohort
   contains 52 mammalian species. It includes
   normalized weights, pre-normalization opportunity-adjusted burdens, and two
   complementary-substitution ratios calculated within each species.
5. `5context192` asks whether the flanking context of four transition types
    changes across the 52-species mammalian three-gene cohort. It keeps structural
   opportunity zeros separate from possible but unobserved contexts. Its
   `signatures` substage also fits COSMIC v3.3 SBS96 activities to strand-aware
   high, low, and high-minus-low transformations of each matched mean gene
   spectrum using both SigProfilerAssignment (Python) and mSigAct (R).
6. `6amino_acid_shift` discovers Mammalia directly in MIDORI2 GB272 and obtains
   exact CDS/protein pairs for ten heavy-strand Major Arc genes, from COX1 to
   CytB, without using a MutSpec species list. Its strict ten-gene QC cohort
   contains 1,932 taxa. The notebook compares `(Asn+Lys)/Gly` and
   `Pro/(Phe+Leu[TTA/TTG])`, draws paired ten-gene boxplots, and displays all 45
   forward gene-pair effects in compact matrices. It does not plot species
   counts or availability.

The main input for stages 2--4 is
`1init_data/data/MutSpecVertebrates12.csv.gz`; stage 5 uses the corresponding
`MutSpecVertebrates192.csv.gz` table. Stage 6 is independent
of those inputs and uses only its cached MIDORI2 sequence archives. Run stages
1--4 first, then either or both stage-5 branches as needed; stage 6 can be run
independently. Stage 6 keeps only a source manifest and the final QC sequence
table on disk; its derived analytical tables remain in notebook memory.

See the README inside each stage for its assumptions, commands, outputs, and
interpretation limits.

The positional proxy used in stages 4 and 5 is calculated from fixed human rCRS
coordinates and is not a species-specific measurement. Restricting the project
to mammals removes between-class pooling but does not remove phylogenetic
dependence, uneven sampling among mammalian clades, or coordinate-system effects.
