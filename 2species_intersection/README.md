# 2. Species intersection

## Purpose

This analysis determines how many mammal species have usable mutation spectra
for the same genes. It reports all pairwise and three-gene intersections and,
when possible, four-gene intersections before any spectra are compared.

## Inputs

The default input is
`../1init_data/data/MutSpecVertebrates12.csv.gz`. Availability means that a
species has a complete, valid 12-component spectrum for a gene, not merely that
a sequence or annotation exists.

## Workflow

`count_common_species.py`:

1. selects the requested taxonomic class (`Mammalia` by default);
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

Use `--taxonomic-class`, `--max-combination-size`, `--input`, or `--output-dir`
to change the defaults. The functions in `count_common_species.py` can also be
imported after adding this stage directory to `sys.path`.

## Outputs

The `data/` folder contains:

- `gene_species_counts.csv`: usable species count for each gene;
- `species_intersection_counts.csv`: intersection size, union size, and Jaccard
  index for each gene combination;
- `species_intersection_members.csv`: one row per species in each intersection.

The `figures/` folder contains `pairwise_common_species_counts.png`.

## Caveats

- Intersections are inclusive: a species present in four genes also contributes
  to every relevant pair and triple.
- Matching uses exact gene and species labels from the source table; it does not
  perform taxonomic synonym resolution.
- The member table, rather than only a reported count, should be used to
  reproduce a downstream matched cohort.
