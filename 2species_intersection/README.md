# 2. Species intersection

## Purpose

This analysis determines how many mammalian species have usable mutation
spectra for the same genes. The current project scope is fixed to `Mammalia`.
It reports all pairwise and three-gene intersections and, when possible,
four-gene intersections before any spectra are compared.

With the supplied input, the mammalian analysis gives 104 common species for
`CO1 + Cytb` and 52 for `CO1 + CO3 + Cytb`. Counts are regenerated
from the input and can change if the source table or validation rules change.

## Inputs

The default input is
`../1init_data/data/MutSpecVertebrates12.csv.gz`. Availability means that a
species has a complete, valid 12-component spectrum for a gene, not merely that
a sequence or annotation exists.

## Workflow

`count_common_species.py`:

1. selects the `Mammalia` rows before profile validation;
2. validates identifiers, the 12 substitution categories, spectrum values, and
   profile sums;
3. builds the set of usable species for each gene;
4. calculates inclusive set intersections for combinations of two through four
   genes; and
5. saves both counts and the exact species members for auditing and reuse.

`CommonSpeciesStatistics.ipynb` presents the same results compactly and plots
the pairwise intersection sizes.

Run the default analysis from the repository root with:

```bash
python 2species_intersection/count_common_species.py
```

Use `--max-combination-size`, `--input`, or `--output-dir` to change the other
defaults. The current taxonomic scope is explicit in the command below:

```bash
python 2species_intersection/count_common_species.py --taxonomic-class Mammalia
```

Omitting `--taxonomic-class` still uses its `Mammalia` default. Matching the
class value is case-insensitive, and no taxonomic aliases are introduced. The
functions in `count_common_species.py` can also be imported after adding this
stage directory to `sys.path`.

## Outputs

The `data/` folder contains:

- `gene_species_counts.csv`: usable species count for each gene;
- `species_intersection_counts.csv`: intersection size, union size, and Jaccard
  index for each gene combination;
- `species_intersection_members.csv`: one row per species in each intersection,
  including its source taxonomic class.

All tables record `Mammalia` as the selected taxonomic scope.

The `figures/` folder contains `pairwise_common_species_counts.png`.

## Verification

`python verify_cohorts.py`, in the repository root, recomputes every saved
intersection from the source table and also checks that the matched cohorts of
stages 3-5 are exactly those intersections.

## Caveats

- Intersections are inclusive: a species present in four genes also contributes
  to every relevant pair and triple.
- Matching uses exact gene and species labels from the source table; it does not
  perform taxonomic synonym resolution.
- Species are not phylogenetically independent, so the counts describe data
  availability rather than independent evolutionary replicates.
- The member table, rather than only a reported count, should be used to
  reproduce a downstream matched cohort.
