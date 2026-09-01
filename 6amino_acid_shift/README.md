# 6. MIDORI2 amino-acid ratios along the Major Arc

## Purpose

This stage compares two amino-acid count ratios across ten mammalian
mitochondrial protein-coding genes in Major Arc order:

`COX1 → COX2 → ATP8 → ATP6 → COX3 → ND3 → ND4L → ND4 → ND5 → CytB`.

The ratios are:

- `(Asn + Lys) / Gly`;
- `Pro / (Phe + Leu[TTA/TTG])`, where only leucines encoded by `TTA` or `TTG`
  enter the denominator.

The less specific `Pro/Phe` ratio is no longer calculated. ND1 and ND2 are
outside the Major Arc definition used by the project. ND6 lies physically in
the arc but is excluded because it is encoded on the opposite, light strand.
MIDORI names ATP6 and ATP8 `A6` and `A8`; the notebook uses their familiar
display names.

Species are discovered directly in MIDORI2. There is no MutSpec intersection,
mutation-spectrum weighting, or GenBank retrieval code.

## MIDORI2 input and cohort

The cache uses MIDORI2 `GenBank272_2026-06-07` (`GB272`) and contains one
`UNIQ_NUC` plus one `TOTAL_AA` archive for each of the ten genes: 20 archives in
total. The large archives live in `data/midori/`, are excluded from Git, and are
not replaced automatically. Their official URLs, sizes, SHA-256 hashes, and
scan counts are stored in `data/derived_midori/midori_file_manifest.csv`.

The backend scans only records annotated as `Mammalia` (TaxID 40674). TaxID is
the taxon key; species names are labels. A protein is paired with a CDS only
when accession and coordinates match exactly. Candidate selection is
deterministic and is completed before any ratio is calculated.

A selected record must be full-length, in frame, nearly unambiguous, free of
internal stops, and reproduce its paired protein under vertebrate mitochondrial
translation table 2. The final analytical cohort requires a QC-passing record
for every target gene, retaining a strictly paired comparison:

- 1,932 mammalian species;
- 19,320 `species × gene` sequence records.

The MIDORI `UNIQ_NUC` representatives are not a strict RefSeq-only collection.

## Counts, zeros, and ratios

P, F, N, K, and G are counted in the validated protein. `Leu[TTA/TTG]` is
counted from the matched CDS because synonymous codons cannot be reconstructed
from a protein sequence. Dividing numerator and denominator by the same gene
length would cancel, so no additional length normalization is used.

Short genes make zero counts important: 1,819 of 1,932 ATP8 records contain no
glycine, and seven ND4L records contain no proline. Raw counts and raw ratios
where defined remain available in notebook memory. Tests, fold changes, and
figures use the explicitly documented continuity correction

```text
corrected ratio = (numerator + 0.5) / (denominator + 0.5)
```

This keeps zero-count observations instead of silently deleting most ATP8
measurements. No analytical CSV tables are generated.

## Statistical comparison

All ten gene measurements are linked within species. For each ratio the
notebook therefore uses:

1. a Friedman repeated-measures test for any difference among the ten genes;
2. Kendall's `W` as the omnibus effect size;
3. a one-sided Page trend test for an overall increase along the predefined
   Major Arc order;
4. all 45 forward pairwise comparisons, tested as `later > earlier` with
   one-sided paired Wilcoxon signed-rank tests on log2 corrected ratios;
5. Holm correction across the 45 post-hoc tests within each ratio.

The nine adjacent pairs and the `CytB / COX1` endpoint are marked as planned,
but all 45 comparisons are calculated. Nothing is selected post hoc because it
looks especially strong. ATP8/ATP6 and ND4L/ND4 overlap in physical coordinates,
so the ordered test describes a broad positional trend rather than ten equally
spaced independent exposure steps.

## Figures

The notebook produces two figures as PNG and editable PDF:

1. `figures/amino_acid_ratios_by_gene.*` — two vertically aligned
   violin/boxplot panels across all ten genes. Faint lines retain every
   within-species trajectory; dots and the dark line show medians. The y-axis is
   logarithmic. As requested, each panel contains only the one-sided paired
   Wilcoxon p-value for the pre-specified `CytB > COX1` endpoint.
2. `figures/paired_ratio_fold_changes.*` — two upper-triangular matrices showing
   the median within-species fold change for every one of the 45 forward gene
   pairs. Adjacent comparisons are outlined, and `CytB / COX1` has a gold
   outline. This avoids both an unreadable 90-row forest plot and post-hoc
   selection of only the largest effects. P-values are not drawn in this figure.

There are no plots of species counts, availability, or attrition.

## Current result

Both ratios differ strongly among genes (`Friedman p < 10^-300`), with
Kendall's `W = 0.927` for `(Asn+Lys)/Gly` and `W = 0.830` for the
proline-based ratio.

The broad ordered Page trend is supported for `(Asn+Lys)/Gly`
(`Holm p = 2.13 × 10^-87`) but not for `Pro/(Phe+Leu[TTA/TTG])` (`p = 1`).
Neither ratio follows a strict step-by-step increase because several adjacent
contrasts decrease. Nevertheless, the pre-specified CytB endpoint is above
COX1 for both ratios:

| Ratio | Median CytB/COX1 fold change | Species with increase |
|---|---:|---:|
| `(Asn+Lys)/Gly` | 1.863× | 100.0% |
| `Pro/(Phe+Leu[TTA/TTG])` | 1.249× | 95.4% |

These are sample-level paired results, not phylogenetically corrected inference.

## Minimal outputs

The extraction backend writes only:

- `data/derived_midori/midori_file_manifest.csv`;
- `data/derived_midori/midori_gene_sequences.csv.gz`.

Counts, ratios, omnibus tests, all 90 ratio-specific post-hoc rows, and effect
matrices exist as in-memory DataFrames inside `AminoAcidShift.ipynb`.

## Execution and verification

Refresh extraction and QC after changing the cached archives:

```powershell
py -3 6amino_acid_shift/run_analysis.py
```

Run `AminoAcidShift.ipynb` from the repository root or stage directory, then:

```powershell
py -3 6amino_acid_shift/verify_midori.py
```

## Interpretation limit

The sample is uneven across mammalian clades, and species are not
phylogenetically independent. Effect sizes and consistency across paired species
are therefore emphasized over extremely small p-values. A phylogenetic or
genus-balanced sensitivity analysis would be a separate extension.
