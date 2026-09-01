"""Extract paired mammalian major-arc H-strand genes from MIDORI2.

This module is intentionally limited to data acquisition, deterministic record
selection, and sequence QC. Statistical analyses and figures live in
``AminoAcidShift.ipynb``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import gzip
import hashlib
from pathlib import Path
import re
from typing import Callable, Iterable, Iterator, Mapping

import pandas as pd

import amino_acid_shift as core


MIDORI_RELEASE = "GenBank272_2026-06-07"
MIDORI_VERSION = "GB272"
MIDORI_BASE_URL = (
    "https://www.reference-midori.info/download/Databases/"
    f"{MIDORI_RELEASE}"
)
MIDORI_FILES: Mapping[str, Mapping[str, str]] = {
    "CO1": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_CO1_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_CO1_RAW_AA.fasta.gz",
    },
    "CO2": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_CO2_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_CO2_RAW_AA.fasta.gz",
    },
    "A8": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_A8_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_A8_RAW_AA.fasta.gz",
    },
    "A6": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_A6_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_A6_RAW_AA.fasta.gz",
    },
    "CO3": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_CO3_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_CO3_RAW_AA.fasta.gz",
    },
    "ND3": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_ND3_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_ND3_RAW_AA.fasta.gz",
    },
    "ND4L": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_ND4L_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_ND4L_RAW_AA.fasta.gz",
    },
    "ND4": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_ND4_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_ND4_RAW_AA.fasta.gz",
    },
    "ND5": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_ND5_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_ND5_RAW_AA.fasta.gz",
    },
    "Cytb": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_Cytb_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_Cytb_RAW_AA.fasta.gz",
    },
}

_CLASS_TAXIDS = {"Mammalia": "40674"}
_SPECIES_TAXON = re.compile(r"(?:^|;)species_([^;]+)_([0-9]+)(?:;|$)")
_CLASS_TAXON = re.compile(r"(?:^|;)class_([^;]+)_([0-9]+)(?:;|$)")
_LOCUS_IDENTIFIER = re.compile(r"^(.+)\.([<>]?\d+)\.([<>]?\d+)$")
_AA_LOCUS_PREFIX = re.compile(r"^(.+\.[<>]?\d+\.[<>]?\d+)_")
_SAFE_SPECIES = re.compile(r"[^A-Za-z0-9_.-]+")


def _normal_name(value: object) -> str:
    return " ".join(str(value).split()).casefold()


def _iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    identifier = ""
    taxonomy = ""
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if identifier:
                    yield identifier, taxonomy, "".join(chunks).upper()
                parts = line[1:].split(maxsplit=1)
                identifier = parts[0]
                taxonomy = parts[1] if len(parts) == 2 else ""
                chunks = []
            elif line:
                chunks.append(line.strip())
    if identifier:
        yield identifier, taxonomy, "".join(chunks).upper()


def _taxonomy_fields(taxonomy: str) -> dict[str, str]:
    species = _SPECIES_TAXON.search(taxonomy)
    class_taxon = _CLASS_TAXON.search(taxonomy)
    return {
        "organism": species.group(1) if species else "",
        "taxid": species.group(2) if species else "",
        "class_name": class_taxon.group(1) if class_taxon else "",
        "class_taxid": class_taxon.group(2) if class_taxon else "",
    }


def _is_requested_class(fields: Mapping[str, str], class_filter: str) -> bool:
    return (
        fields["class_name"] == class_filter
        or fields["class_taxid"] == _CLASS_TAXIDS[class_filter]
    )


def _parse_locus_identifier(identifier: str) -> dict[str, object]:
    match = _LOCUS_IDENTIFIER.fullmatch(identifier)
    if not match:
        return {
            "accession_version": "",
            "accession": "",
            "start_text": "",
            "end_text": "",
            "is_partial": True,
        }
    accession_version, start_text, end_text = match.groups()
    return {
        "accession_version": accession_version,
        "accession": re.sub(r"\.\d+$", "", accession_version),
        "start_text": start_text,
        "end_text": end_text,
        "is_partial": start_text.startswith("<") or end_text.startswith(">"),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def midori_file_manifest(
    cache_dir: str | Path,
    scan_counts: Mapping[tuple[str, str], Mapping[str, int]] | None = None,
) -> pd.DataFrame:
    cache = Path(cache_dir)
    rows = []
    for gene in core.TARGET_GENES:
        for molecule in ("NUC", "AA"):
            relative_url = MIDORI_FILES[gene][molecule]
            path = cache / Path(relative_url).name
            if not path.exists():
                raise FileNotFoundError(
                    f"Missing MIDORI2 input {path}. Download it from "
                    f"{MIDORI_BASE_URL}/{relative_url}"
                )
            rows.append(
                {
                    "release": MIDORI_RELEASE,
                    "version": MIDORI_VERSION,
                    "database_type": "UNIQ" if molecule == "NUC" else "TOTAL",
                    "Gene": gene,
                    "molecule": molecule,
                    "filename": path.name,
                    "url": f"{MIDORI_BASE_URL}/{relative_url}",
                    "compressed_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                    **dict((scan_counts or {}).get((gene, molecule), {})),
                }
            )
    return pd.DataFrame(rows)


def _scan_nucleotides(
    cache: Path,
    class_filter: str,
    progress: Callable[[str], None] | None,
) -> tuple[
    dict[str, dict[str, list[dict[str, object]]]],
    dict[str, Counter[str]],
    dict[str, set[str]],
    dict[tuple[str, str], dict[str, int]],
]:
    by_gene: dict[str, dict[str, list[dict[str, object]]]] = {}
    names_by_taxid: dict[str, Counter[str]] = defaultdict(Counter)
    genes_by_taxid: dict[str, set[str]] = defaultdict(set)
    scan_counts: dict[tuple[str, str], dict[str, int]] = {}

    for gene in core.TARGET_GENES:
        path = cache / Path(MIDORI_FILES[gene]["NUC"]).name
        if not path.exists():
            midori_file_manifest(cache)
        if progress:
            progress(f"Scanning {path.name}")
        records: dict[str, list[dict[str, object]]] = defaultdict(list)
        n_records = 0
        n_class = 0
        for identifier, taxonomy, sequence in _iter_fasta(path):
            n_records += 1
            fields = _taxonomy_fields(taxonomy)
            if not _is_requested_class(fields, class_filter):
                continue
            n_class += 1
            taxid = fields["taxid"]
            organism = fields["organism"]
            if not taxid or not organism:
                continue
            records[taxid].append(
                {
                    "identifier": identifier,
                    "sequence": sequence,
                    **fields,
                }
            )
            names_by_taxid[taxid][organism] += 1
            genes_by_taxid[taxid].add(gene)
        by_gene[gene] = dict(records)
        scan_counts[(gene, "NUC")] = {
            "n_records": n_records,
            "n_requested_class": n_class,
            "n_retained_candidates": sum(map(len, records.values())),
            "n_retained_taxa": len(records),
        }
    return by_gene, names_by_taxid, genes_by_taxid, scan_counts


def _scan_proteins(
    cache: Path,
    class_filter: str,
    nucleotides: Mapping[str, Mapping[str, list[dict[str, object]]]],
    scan_counts: dict[tuple[str, str], dict[str, int]],
    progress: Callable[[str], None] | None,
) -> dict[str, dict[str, list[dict[str, object]]]]:
    protein_index: dict[str, dict[str, list[dict[str, object]]]] = {}
    for gene in core.TARGET_GENES:
        wanted_loci = {
            str(record["identifier"]): taxid
            for taxid, records in nucleotides[gene].items()
            for record in records
        }
        path = cache / Path(MIDORI_FILES[gene]["AA"]).name
        if not path.exists():
            midori_file_manifest(cache)
        if progress:
            progress(f"Scanning {path.name}")
        records: dict[str, list[dict[str, object]]] = defaultdict(list)
        retained_taxids: set[str] = set()
        n_records = 0
        n_class = 0
        for identifier, taxonomy, sequence in _iter_fasta(path):
            n_records += 1
            fields = _taxonomy_fields(taxonomy)
            if not _is_requested_class(fields, class_filter):
                continue
            n_class += 1
            match = _AA_LOCUS_PREFIX.match(identifier)
            if not match:
                continue
            locus_identifier = match.group(1)
            expected_taxid = wanted_loci.get(locus_identifier)
            if expected_taxid is None or fields["taxid"] != expected_taxid:
                continue
            records[locus_identifier].append(
                {"identifier": identifier, "sequence": sequence, **fields}
            )
            retained_taxids.add(fields["taxid"])
        protein_index[gene] = dict(records)
        scan_counts[(gene, "AA")] = {
            "n_records": n_records,
            "n_requested_class": n_class,
            "n_retained_candidates": sum(map(len, records.values())),
            "n_retained_taxa": len(retained_taxids),
        }
    return protein_index


def _analysis_catalog(
    names_by_taxid: Mapping[str, Counter[str]],
    complete_taxids: set[str],
) -> pd.DataFrame:
    canonical_names = {
        taxid: sorted(names_by_taxid[taxid].items(), key=lambda item: (-item[1], item[0]))[0][0]
        for taxid in complete_taxids
    }
    slugs = {
        taxid: _SAFE_SPECIES.sub("_", name.replace(" ", "_")).strip("_")
        for taxid, name in canonical_names.items()
    }
    duplicate_slugs = Counter(slugs.values())
    rows = []
    for taxid in sorted(complete_taxids, key=lambda value: (canonical_names[value], value)):
        slug = slugs[taxid]
        species = slug if duplicate_slugs[slug] == 1 else f"{slug}__taxid_{taxid}"
        rows.append(
            {
                "Species": species,
                "Species_query": canonical_names[taxid],
                "TaxID": taxid,
            }
        )
    return pd.DataFrame(rows)


def _matching_protein(
    nucleotide: Mapping[str, object],
    proteins: Iterable[Mapping[str, object]],
) -> tuple[dict[str, object] | None, str]:
    prefix = f"{nucleotide['identifier']}_"
    paired = {
        str(record["identifier"]): dict(record)
        for record in proteins
        if str(record["identifier"]).startswith(prefix)
        and str(record.get("taxid", "")) == str(nucleotide.get("taxid", ""))
    }
    if len(paired) == 1:
        return next(iter(paired.values())), "same_accession_and_coordinates"
    if len(paired) > 1:
        return None, "ambiguous_same_locus_AA"
    return None, "no_same_locus_AA"


def _translation_parameters(cds: str, protein: str) -> tuple[int, bool, bool]:
    protein = protein.rstrip("*")
    for codon_start in (1, 2, 3):
        coding = cds[codon_start - 1:]
        for apply_start_rule in (True, False):
            calculated = core.translate_table2(
                coding,
                apply_start_rule=apply_start_rule,
            ).rstrip("*")
            if protein and calculated == protein:
                return codon_start, apply_start_rule, True
    return 1, True, False


def _prepare_candidate(
    target: Mapping[str, object],
    gene: str,
    nucleotide: Mapping[str, object],
    proteins: Iterable[Mapping[str, object]],
) -> dict[str, object]:
    protein, pair_status = _matching_protein(nucleotide, proteins)
    locus = _parse_locus_identifier(str(nucleotide["identifier"]))
    cds = str(nucleotide["sequence"])
    protein_sequence = str(protein["sequence"]) if protein else ""
    codon_start, apply_start_rule, translation_match = _translation_parameters(
        cds,
        protein_sequence,
    )
    coding = cds[codon_start - 1:]
    remainder = len(coding) % 3
    terminal_fragment = coding[-remainder:] if remainder else ""
    frame_valid = remainder == 0 or terminal_fragment in {"T", "TA"}
    calculated = core.translate_table2(coding, apply_start_rule=apply_start_rule)
    ambiguous = sum(base not in "ACGT" for base in coding)
    ambiguous_fraction = ambiguous / len(coding) if coding else 1.0
    protein_length = len(protein_sequence.rstrip("*"))
    low, high = core.PROTEIN_LENGTH_BOUNDS[gene]
    basic_valid = bool(
        protein
        and pair_status == "same_accession_and_coordinates"
        and translation_match
        and not bool(locus["is_partial"])
        and ambiguous_fraction <= 0.01
        and frame_valid
        and "*" not in calculated[:-1]
        and low <= protein_length <= high
    )
    is_refseq = str(locus["accession_version"]).startswith(("NC_", "NW_", "NZ_"))
    sort_key = (
        0 if basic_valid else 1,
        0 if is_refseq else 1,
        0 if pair_status == "same_accession_and_coordinates" else 1,
        0 if translation_match else 1,
        0 if not bool(locus["is_partial"]) else 1,
        ambiguous_fraction,
        abs(protein_length - (low + high) / 2),
        str(nucleotide["identifier"]),
    )
    row = {
        "Species": target["Species"],
        "Species_query": target["Species_query"],
        "TaxID": target["TaxID"],
        "Gene": gene,
        "accession": locus["accession"],
        "accession_version": locus["accession_version"],
        "MIDORI_organism": nucleotide["organism"],
        "coding_strand": 1,
        "transl_table": 2,
        "cds_sequence": cds,
        "protein_sequence": protein_sequence,
        "cds_length": len(cds),
        "protein_length": protein_length,
        "cds_coordinates": f"{locus['start_text']}..{locus['end_text']}",
        "codon_start": codon_start,
        "apply_start_rule": apply_start_rule,
        "is_partial": locus["is_partial"],
        "ambiguous_nucleotides": ambiguous,
        "is_refseq_accession": is_refseq,
        "midori_nuc_identifier": nucleotide["identifier"],
        "midori_aa_identifier": protein["identifier"] if protein else "",
        "dna_aa_pair_status": pair_status,
        "same_feature_pair": pair_status == "same_accession_and_coordinates",
        "translation_match_inferred": translation_match,
        "midori_release": MIDORI_RELEASE,
        "selection_reason": "",
    }
    return {
        "accession_version": str(locus["accession_version"]),
        "basic_valid": basic_valid,
        "sort_key": sort_key,
        "row": row,
    }


def retrieve_midori_sequences(
    cache_dir: str | Path,
    *,
    class_filter: str = "Mammalia",
    progress: Callable[[str], None] | None = print,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Return one selected record per complete MIDORI taxon and target gene."""

    if class_filter not in _CLASS_TAXIDS:
        raise ValueError("Only Mammalia is currently supported")
    cache = Path(cache_dir)
    nucleotides, names, genes_by_taxid, scan_counts = _scan_nucleotides(
        cache,
        class_filter,
        progress,
    )
    complete_taxids = {
        taxid for taxid, genes in genes_by_taxid.items()
        if set(core.TARGET_GENES).issubset(genes)
    }
    nucleotides = {
        gene: {
            taxid: records for taxid, records in gene_records.items()
            if taxid in complete_taxids
        }
        for gene, gene_records in nucleotides.items()
    }
    proteins = _scan_proteins(
        cache,
        class_filter,
        nucleotides,
        scan_counts,
        progress,
    )
    catalog = _analysis_catalog(names, complete_taxids)

    selected_rows: list[dict[str, object]] = []
    for target in catalog.to_dict("records"):
        taxid = str(target["TaxID"])
        candidates_by_gene = {
            gene: sorted(
                [
                    _prepare_candidate(
                        target,
                        gene,
                        nucleotide,
                        proteins[gene].get(str(nucleotide["identifier"]), ()),
                    )
                    for nucleotide in nucleotides[gene][taxid]
                ],
                key=lambda candidate: candidate["sort_key"],
            )
            for gene in core.TARGET_GENES
        }
        valid_accessions = {
            gene: {
                candidate["accession_version"]
                for candidate in candidates
                if candidate["basic_valid"] and candidate["accession_version"]
            }
            for gene, candidates in candidates_by_gene.items()
        }
        shared_accessions = set.intersection(
            *(valid_accessions[gene] for gene in core.TARGET_GENES)
        )
        if shared_accessions:
            shared_accession = min(
                shared_accessions,
                key=lambda accession: (
                    0 if accession.startswith(("NC_", "NW_", "NZ_")) else 1,
                    accession,
                ),
            )
            selected = {
                gene: next(
                    candidate for candidate in candidates_by_gene[gene]
                    if candidate["basic_valid"]
                    and candidate["accession_version"] == shared_accession
                )
                for gene in core.TARGET_GENES
            }
            reason = "valid exact DNA-AA pairs from one shared accession"
        else:
            selected = {
                gene: candidates_by_gene[gene][0]
                for gene in core.TARGET_GENES
            }
            reason = "best deterministic gene-specific candidate"
        accessions = [selected[gene]["accession_version"] for gene in core.TARGET_GENES]
        same_accession = len(set(accessions)) == 1
        for gene in core.TARGET_GENES:
            row = dict(selected[gene]["row"])
            row["same_accession_all_targets"] = same_accession
            row["selection_reason"] = reason
            selected_rows.append(row)

    summary = {
        "n_taxa_any_target_gene": len(names),
        "n_taxa_all_target_nucleotide_genes": len(complete_taxids),
    }
    return pd.DataFrame(selected_rows), midori_file_manifest(cache, scan_counts), summary


def run_midori_analysis(
    stage_dir: str | Path,
    *,
    class_filter: str = "Mammalia",
    progress: Callable[[str], None] | None = print,
) -> dict[str, object]:
    """Extract MIDORI data, apply QC, and persist only the essential inputs."""

    stage = Path(stage_dir)
    cache = stage / "data" / "midori"
    derived = stage / "data" / "derived_midori"
    cache.mkdir(parents=True, exist_ok=True)
    derived.mkdir(parents=True, exist_ok=True)

    selected, file_manifest, summary = retrieve_midori_sequences(
        cache,
        class_filter=class_filter,
        progress=progress,
    )
    qc = core.sequence_qc(selected)
    n_target_genes = len(core.TARGET_GENES)
    exact_pair_species = qc.groupby("Species")["same_feature_pair"].agg(
        lambda values: len(values) == n_target_genes and bool(values.all())
    )
    valid_species = qc.groupby("Species")["sequence_qc_status"].agg(
        lambda values: len(values) == n_target_genes and bool(values.eq("pass").all())
    )
    analysis_species = set(valid_species.index[valid_species])
    analysis_qc = qc.loc[qc["Species"].isin(analysis_species)].copy()
    analysis_qc = analysis_qc.sort_values(["Species", "Gene"], kind="stable").reset_index(drop=True)

    summary.update(
        {
            "n_taxa_all_target_exact_pairs": int(exact_pair_species.sum()),
            "n_taxa_final_paired_qc": len(analysis_species),
        }
    )
    # Kept in memory for backwards compatibility with early copies of the
    # notebook, which display ``results['attrition']``.  This compact view is
    # not written as another derived table and does not restore the removed
    # availability/attrition figures.
    attrition = pd.DataFrame(
        {
            "selection_step": [
                "at least one target gene",
                "all target nucleotide genes",
                "exact CDS/protein pairs for all targets",
                "final paired sequence-QC cohort",
            ],
            "n_taxa": [
                summary["n_taxa_any_target_gene"],
                summary["n_taxa_all_target_nucleotide_genes"],
                summary["n_taxa_all_target_exact_pairs"],
                summary["n_taxa_final_paired_qc"],
            ],
        }
    )
    attrition["retained_from_previous"] = (
        attrition["n_taxa"] / attrition["n_taxa"].shift(1)
    )

    file_manifest.to_csv(derived / "midori_file_manifest.csv", index=False)
    analysis_qc.to_csv(
        derived / "midori_gene_sequences.csv.gz",
        index=False,
        compression="gzip",
    )
    return {
        "file_manifest": file_manifest,
        "summary": summary,
        "attrition": attrition,
        "analysis_qc": analysis_qc,
        "retrieval_backend": f"MIDORI2 {MIDORI_VERSION} UNIQ_NUC+TOTAL_AA RAW",
        "taxonomic_scope": class_filter,
    }


__all__ = [
    "MIDORI_BASE_URL",
    "MIDORI_FILES",
    "MIDORI_RELEASE",
    "MIDORI_VERSION",
    "midori_file_manifest",
    "retrieve_midori_sequences",
    "run_midori_analysis",
]
