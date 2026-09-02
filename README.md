# TBSS-Signature

Does the mitochondrial mutation spectrum, and the protein composition it
leaves behind, change along a positional proxy for the time that mtDNA spends
single-stranded during replication?

The taxonomic scope is `Mammalia` throughout. The supplied mutation tables keep
other vertebrate classes as immutable source data, but every analysis filters to
mammals before forming a cohort. Each numbered folder is a self-contained stage
with its own notebook or script, generated tables, figures, and a README that
documents that stage's assumptions and limits.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv/Scripts/activate
python -m pip install -r requirements.txt
python -m pip install jupyterlab
```

Run the stages in order; every notebook works from the repository root or from
its own stage folder. Then check the two data contracts:

```bash
python verify_cohorts.py                    # stages 2-5: cohorts and gene order
python 6amino_acid_shift/verify_midori.py   # stage 6: MIDORI2 data and figures
```

Both verifiers are plain scripts with no test framework and no network access.
The COSMIC signature substage in `5context192/signatures/` needs its own
environment and R; see its README.

## Two data sources, two species sets

The project has two independent inputs, and they are **not** expected to cover
the same species:

| | Stages 1-5 | Stage 6 |
|---|---|---|
| Input | `MutSpecVertebrates12/192.csv.gz` | MIDORI2 GB272 FASTA archives |
| Quantity | reconstructed substitutions from a phylogeny | amino-acid and codon counts |
| Species | 828 mammals, 4 genes (`CO1`, `CO3`, `Cytb`, `ND2`) | 1,932 mammals, 10 genes |

Stages 1-5 can only use species for which mutations were reconstructed on a
tree, whereas stage 6 uses every mammalian record with a usable CDS/protein
pair. The two cohorts therefore differ by construction; this is accepted at the
current stage of the project and no analysis mixes them.

**Within** each input the rule is strict: a gene comparison uses only species
that have a complete, valid profile for *every* selected gene, so a change in
species composition can never masquerade as a gene effect. `verify_cohorts.py`
re-derives each saved cohort from the source table and fails if it is not the
exact intersection.

## Canonical gene order

[`mtdna.py`](mtdna.py) is the single source of mitochondrial gene order,
human rCRS coordinates, strand, display names, and the DssH positional proxy.
Mammalian gene order is conserved, so the repository fixes one order and derives
it nowhere else:

```text
ND1 ND2 | CO1 CO2 A8 A6 CO3 ND3 ND4L ND4 ND5 ND6 | Cytb
```

`canonical_order()` reorders any gene selection into rCRS order, so a figure or
an ordered trend test can never depend on how a gene list was typed in a
notebook, and it rejects unknown gene keys. Stage 6's ten heavy-strand Major Arc
genes and the DssH metadata of stages 4-5 are both derived from this table.

The proxy is `2 * (gene midpoint - OriL) / 16569` with `OriL = 5730`. It is a
linear, unwrapped coordinate: it is meaningful only downstream of OriL, it
exceeds one for CytB, and it is negative for ND1 and ND2, which is why those two
genes cannot be placed on the current gradient without redefining it.

## Repository structure

```text
mtdna.py                          canonical gene order, coordinates, DssH proxy
verify_cohorts.py                 matched-cohort and gene-order contract check
requirements.txt                  pinned Python environment for stages 1-6
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
|-- data/midori/                  cached MIDORI2 GB272 FASTA archives (not in git)
|-- data/derived_midori/          source manifest and final QC sequences
|-- figures/                      ratio distributions and paired effects
|-- midori_analysis.py            Mammalia discovery, pairing, and selection
|-- amino_acid_shift.py           translation-table and sequence-QC helpers
|-- utils.py                      shared constants and checksums
|-- run_analysis.py               MIDORI extraction/QC entry point
|-- verify_midori.py              stage-6 data and figure contract check
|-- AminoAcidShift.ipynb
`-- README.md
```

## Analysis order

1. `1init_data` documents and checks the supplied mutation-spectrum tables.
2. `2species_intersection` counts mammalian species shared by every gene
   combination and records the exact members of each matched cohort.
3. `3compare_t_genes` compares 12-component spectra for two to five genes using
   the same species for every selected gene. The current `CO1`/`Cytb` example
   contains 104 mammalian species and includes paired boxplots for `C>T`, `G>A`,
   and the `C>T/G>A` and `A>G/T>C` within-profile ratios, with pre-specified
   one-sided paired Wilcoxon tests.
4. `4tsss_gradient` explores changes in selected mutation types along the DssH
   proxy. Its three-gene cohort contains 52 mammalian species. It includes
   normalized weights, pre-normalization opportunity-adjusted burdens, and two
   complementary-substitution ratios calculated within each species.
5. `5context192` asks whether the flanking context of four transition types
   changes across the same 52-species three-gene cohort. It keeps structural
   opportunity zeros separate from possible but unobserved contexts. Its
   `signatures` substage fits COSMIC v3.3 SBS96 activities to strand-aware
   high, low, and high-minus-low transformations of each matched mean gene
   spectrum using both SigProfilerAssignment (Python) and mSigAct (R).
6. `6amino_acid_shift` discovers Mammalia directly in MIDORI2 GB272 and obtains
   exact CDS/protein pairs for the ten heavy-strand Major Arc genes, from COX1
   to CytB, without using a MutSpec species list. Its strict ten-gene QC cohort
   contains 1,932 taxa. The notebook compares `(Asn+Lys)/Gly` and
   `Pro/(Phe+Leu[TTA/TTG])`, draws paired ten-gene boxplots, displays all 45
   forward gene-pair effects in compact matrices, and repeats every test
   separately without the 0.5 continuity correction on the genes that never
   need it.

The main input for stages 2-4 is `1init_data/data/MutSpecVertebrates12.csv.gz`;
stage 5 uses `MutSpecVertebrates192.csv.gz`. Stage 6 is independent of those
inputs and uses only its cached MIDORI2 archives, so it can be run on its own.

## What is tracked

Source code, source data, the derived tables each README refers to, and every
published figure are in git. Downloaded archives (the ~700 MB MIDORI2 FASTA
cache), virtual environments, R package libraries, and scratch files that a
script rebuilds are not; see the comments in [`.gitignore`](.gitignore).
`.gitattributes` stores all text with LF so Windows and Linux checkouts produce
identical diffs.

Stage 6 keeps only a source manifest and the final QC sequence table on disk;
its derived analytical tables stay in notebook memory by design.

## Interpretation limits

The positional proxy is calculated from fixed human rCRS coordinates and is not
a species-specific measurement. Restricting the project to mammals removes
between-class pooling but does not remove phylogenetic dependence, uneven
sampling among mammalian clades, or coordinate-system effects. Species are
treated as independent observations in every test, so p-values are descriptive
evidence for the sampled cohort rather than phylogenetically corrected
inference.

See the README inside each stage for its assumptions, commands, outputs, and
interpretation limits.
