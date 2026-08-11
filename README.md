# TBSS-Signature

This project studies whether mitochondrial mutation spectra vary along a
positional proxy for the time that mtDNA remains single-stranded. The default
scope is all vertebrate classes represented in the supplied table:
`Actinopteri`, `Amphibia`, `Aves`, `Lepidosauria`, and `Mammalia`. The repository
is organised as a sequence of self-contained analyses. Each numbered folder
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
5amino_acid_shift/
|-- data/                         theoretical and spectrum-weighted AA tables
|-- figures/                      amino-acid shift plots
|-- amino_acid_shift.py           reusable genetic-code and bootstrap functions
|-- AminoAcidShift.ipynb
`-- README.md
5context192/
|-- data/                         matched context spectra and summaries
|-- figures/                      192-component and context-gradient plots
|-- context192_analysis.py        reusable context-spectrum functions
|-- Context192MatchedGenes.ipynb
`-- README.md
6amino_acid_shift/
|-- data/midori/                  cached MIDORI2 GB272 FASTA archives
|-- data/derived_midori/          MIDORI sequence, QC, codon, and AA results
|-- figures/midori/               MIDORI mechanistic figure
|-- midori_analysis.py            MIDORI parsing, pairing, and selection
|-- amino_acid_shift.py           shared QC, codon, bootstrap, and plot helpers
|-- run_analysis.py               command-line entry point
|-- AminoAcidShift.ipynb
`-- README.md
```

## Analysis order

1. `1init_data` documents and checks the supplied mutation-spectrum tables.
2. `2species_intersection` counts the vertebrate species shared by every gene
   combination and records the exact members and source class of each matched
   cohort. A single source `Class` value can be selected optionally.
3. `3compare_t_genes` compares 12-component spectra for two to five genes using
   the same species for every selected gene. The current `CO1`/`Cytb` example
   contains 178 species across all five represented classes.
4. `4tsss_gradient` explores changes in selected mutation types along an
   explicit proxy for single-stranded duration. Its current three-gene cohort
   contains 85 species across all five represented classes. It includes
   normalized weights, pre-normalization opportunity-adjusted burdens, and two
   complementary-substitution ratios calculated within each species.
5. Two parallel stage-5 branches extend the three-gene gradient:
   `5amino_acid_shift` uses the 12-component spectra and mitochondrial genetic
   code to estimate theoretical spectrum-weighted amino-acid pressure for the
   exact Mammalia intersection, while `5context192` asks whether the flanking
   context of four transition types changes across all represented classes.
   The amino-acid branch is explicitly theoretical in the absence of CDS/codon
   usage; the context branch keeps structural opportunity zeros separate from
   possible but unobserved contexts and includes pooled and equal-class
   sensitivity summaries.
6. `6amino_acid_shift` obtains real COX1, COX3, and CytB CDS/protein pairs from
   MIDORI2 GB272 for the exact Mammalia mutation-spectrum intersection. Exact
   accession-coordinate identifiers prevent mixing nucleotide and protein
   records, while source/current NCBI TaxIDs verify renamed taxa. Sequence QC,
   observed composition, real-codon opportunities, MutSpec weighting, and
    predicted-versus-observed contrasts use one final matched cohort. The current
    strict MIDORI run retains 51 of 52 source species.

The main input for stages 2--4 and the amino-acid stage-5 branch is
`1init_data/data/MutSpecVertebrates12.csv.gz`; the context stage-5 branch uses
the corresponding `MutSpecVertebrates192.csv.gz` table. Stage 6 also uses the
species taxonomy in `1init_data/data/info.csv` and cached MIDORI2 sequence
archives. Run stages 1--4 first, then either or both stage-5 branches and stage
6 as needed. Tables are written to the corresponding `data/` folder and plots
to `figures/`; source data remain in `1init_data/data/`.

See the README inside each stage for its assumptions, commands, outputs, and
interpretation limits.

The positional proxy used in stages 4 and 5 is calculated from fixed human rCRS
coordinates and is not a species-specific measurement. This limitation is
especially important for the pooled vertebrate analysis: the five classes are
evolutionarily diverse and unevenly represented, so a pooled trend can combine
gene, class, phylogenetic, and coordinate-system effects.
