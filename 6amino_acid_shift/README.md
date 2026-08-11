# 6. MIDORI2 sequences and amino-acid shift

## Purpose

This stage links the exact Mammalia mutation-spectrum intersection for `CO1`,
`CO3`, and `Cytb` to real mitochondrial coding sequences from MIDORI2. It then
uses the matched CDS/protein data for observed amino-acid composition,
real-codon opportunities for `A_H>G_H` and `C_H>T_H`, mutation-spectrum
weighting, and an exploratory predicted-versus-observed comparison.

MIDORI2 is the only sequence backend used by this stage.

## Source release and local cache

The current analysis uses MIDORI2 release `GenBank272_2026-06-07` (`GB272`):

- nucleotide: `RAW/uniq/MIDORI2_UNIQ_NUC_GB272_{gene}_RAW.fasta.gz`;
- protein: `RAW_AA/total/MIDORI2_TOTAL_AA_GB272_{gene}_RAW_AA.fasta.gz`;
- `{gene}` is `CO1`, `CO3`, or `Cytb`.

The six archives are cached under `data/midori/`. Their exact official URL,
compressed size, and SHA-256 checksum are written to
`data/derived_midori/midori_file_manifest.csv`. The analysis never downloads
or replaces these files automatically.

The MIDORI website currently presents an expired TLS certificate. The cached
files were downloaded from the exact paths listed on its official download
page, and local hashes are recorded for repeatability. A local hash detects a
changed cache file but cannot retrospectively authenticate a download made
while TLS verification was unavailable.

Biopython is used here only for resolving merged/current NCBI TaxIDs when the
taxonomy cache is absent. It is not used to retrieve the nucleotide or protein
sequences.

## Why UNIQ_NUC plus TOTAL_AA

MIDORI's `LONGEST_NUC` and `LONGEST_AA` products choose one record per species
independently. In a strict trial, most selected nucleotide loci therefore did
not share an accession-coordinate identifier with the selected protein. In
addition, “longest” records sometimes contained many ambiguous bases.

The implemented route instead scans `UNIQ_NUC` candidates and looks up the
protein in `TOTAL_AA`. A protein is accepted only when its FASTA identifier
starts with the complete nucleotide identifier:

```text
accession.version.start.end_protein_accession.version
```

Thus a protein from the same species but a different nucleotide locus is never
substituted. `translate(CDS)` is checked against the paired MIDORI protein
under vertebrate mitochondrial translation table 2.

## Taxonomy and deterministic selection

Repository species identifiers are converted from underscores to spaces.
Only whitespace and case are normalized. Subspecies, strain, `cf.`, `aff.`,
and `sp.` tokens are not removed. A record must match the exact scientific
name or an NCBI-resolved TaxID; ambiguous/fuzzy matching is not used.

For each species, all MIDORI2 candidates are retained in
`midori_sequence_manifest.csv`. Selection follows:

1. require an exact nucleotide-locus/protein pair;
2. require exact translation, a non-partial CDS, at most 1% ambiguous bases,
   no internal stop, a valid frame or terminal `T`/`TA` incomplete stop, and a
   plausible gene-specific protein length;
3. prefer one accession that supplies basic-valid CO1, CO3, and Cytb pairs;
4. prefer a RefSeq-style accession, fewer ambiguous bases, then a lexical
   identifier as deterministic tie-breaks;
5. if no valid shared accession exists, select the best candidate per gene and
   retain the fallback status in the audit tables.

`UNIQ_NUC` collapses identical nucleotide haplotypes, so a shared accession can
only be selected when that accession remains represented. The analysis reports
`same_accession_all_three` rather than assuming it.

MIDORI FASTA sequences are already emitted in coding orientation. The stored
`coding_strand=1` refers to that supplied sequence orientation; the original
genomic feature strand is not present in MIDORI RAW FASTA.

## Execution

From the repository root:

```powershell
py -3 6amino_acid_shift/run_analysis.py
```

Use `--bootstrap 2000` (the default) for final species-cluster bootstrap
intervals. The analysis is intentionally restricted to Mammalia because its QC
and codon-consequence calculations use vertebrate mitochondrial translation
table 2.

## Current cohort result

For the Mammalia mutation-spectrum input:

- mutation spectra for all three genes: 52 species;
- reliable MIDORI2 sequences for all three genes: 52 species;
- all three genes passed strict sequence QC: 51 species;
- 50 species used one shared accession for all three genes;
- `Ochotona_cansus` passed using a documented gene-specific fallback;
- `Cricetulus_kamensis` was excluded because every available CO3 candidate was
  marked partial at the 3′ boundary.

The exact counts are generated programmatically; they are not used as analysis
constants.

## Outputs

MIDORI2 outputs are isolated under `data/derived_midori/`:

- `midori_file_manifest.csv`: release, URL, size, hash, and scan counts;
- `midori_sequence_manifest.csv`: every candidate, QC prerequisites, rank,
  selected flag, and selection reason;
- `species_midori_matching.csv`: exact/TaxID-confirmed species matching;
- `midori_record_consistency.csv`: three accessions and shared-record status;
- `midori_gene_sequences.csv.gz`: canonical selected CDS/protein table;
- `sequence_qc.csv`;
- `species_gene_availability.csv`;
- `species_intersection_summary.csv`;
- `species_exclusion_reasons.csv` and `exclusion_reason_summary.csv`;
- `final_common_species.csv`;
- observed composition, codon-opportunity, mutation-weighted shift, bootstrap,
  contrast, and Spearman tables.

The figure is written to `figures/midori/`:

- `main_amino_acid_shift.png`.

All confidence intervals resample species as matched clusters. The association
between predicted and observed amino-acid changes is exploratory: gene
identity, selection, phylogeny, base composition, and the gene-order/DssH proxy
remain confounded, so it is not evidence of causality.
