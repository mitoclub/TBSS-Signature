# TBSS-Signature

The aim of the project is to find signatures associated with the time mitochondrial DNA remains single-stranded.

## Reusable matched-spectrum plots

`scripts/mutation_comparison.py` contains strict plotting helpers for two to
five genes. Plotting functions expect an input DataFrame that has already been
restricted to the exact same species for every gene:

```python
from scripts.mutation_comparison import (
    orient_substitutions_to_heavy_strand,
    plot_matched_spectra,
    select_common_species,
)

matched, common_species = select_common_species(mammals, ("CO1", "Cytb"))
matched = orient_substitutions_to_heavy_strand(matched)
result = plot_matched_spectra(matched, genes=("CO1", "Cytb"))
```

The plotter validates the shared cohort and all 12 components but never changes
the cohort itself. See the compact examples in:

- `notebooks/CommonSpeciesStatistics.ipynb`
- `notebooks/CompareTGenes.ipynb`
- `notebooks/TSSSMutationGradient.ipynb`
