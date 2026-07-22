# 1. Initial data

## Purpose

This stage stores the supplied datasets and provides a lightweight quality
check before downstream analyses. Files in `data/` are treated as source inputs:
later stages read them but do not modify them.

## Inputs

- `data/MutSpecVertebrates12.csv.gz` is the primary table of context-free,
  12-component mutation spectra used by stages 2--4.
- `data/MutSpecVertebrates192.csv.gz` contains the context-dependent,
  192-component spectra retained from the earlier workflow.
- `data/obs_muts.csv`, `data/expected.csv`, and `data/exp_cxt_freqs.csv` contain
  supporting observed and expected mutation information.
- `data/info.csv` contains supporting sample metadata.

The downstream 12-component workflow expects at least the columns `Gene`,
`Class`, `Species`, `Mut`, and `MutSpec`. A usable gene/species profile must have
one row for every one of the 12 single-base substitutions, finite non-negative
weights, and a spectrum that sums to one.

## Workflow

`Check192spec.ipynb` is a compact exploratory quality-control notebook retained
from the initial analysis. Despite its historical name, its current input is
the 12-component table. It inspects available genes and spectra and provides
quick visual checks.

Open and run the notebook from either the repository root or this folder:

```text
1init_data/Check192spec.ipynb
```

## Outputs

This stage does not define canonical derived tables. Notebook displays are
diagnostic; authoritative species intersections are generated in
`../2species_intersection/`.

## Caveats

- Gene and species identifiers are preserved exactly as supplied; no taxonomy
  aliases or spelling variants are reconciled here.
- Missing values or incomplete profiles must not be interpreted as biological
  absence.
- The source tables should be replaced only deliberately, because doing so can
  change every downstream cohort and result.
