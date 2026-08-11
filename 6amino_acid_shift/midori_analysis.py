"""MIDORI2 sequence backend for the stage-6 amino-acid analysis.

The adapter reads gene-specific MIDORI2 UNIQ nucleotide and TOTAL amino-acid
RAW FASTA files.  A protein is accepted only when its FASTA identifier starts
with the complete nucleotide locus identifier, which keeps the nucleotide CDS
and amino-acid translation tied to the same accession and coordinates.

MIDORI2 UNIQ_NUC can contain multiple records per species and gene. Selection is
deterministic and first prefers a single accession that supplies all three
valid DNA--AA locus pairs; otherwise it chooses the best valid record per gene.
"""

from __future__ import annotations

from collections import defaultdict
import gzip
import hashlib
from pathlib import Path
import re
from typing import Callable, Iterable, Iterator, Mapping

import numpy as np
import pandas as pd

import amino_acid_shift as core


MIDORI_RELEASE = "GenBank272_2026-06-07"
MIDORI_VERSION = "GB272"
MIDORI_DATABASE_TYPE = "UNIQ_NUC+TOTAL_AA"
MIDORI_BASE_URL = (
    "https://www.reference-midori.info/download/Databases/"
    f"{MIDORI_RELEASE}"
)

MIDORI_FILES: Mapping[str, Mapping[str, str]] = {
    "CO1": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_CO1_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_CO1_RAW_AA.fasta.gz",
    },
    "CO3": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_CO3_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_CO3_RAW_AA.fasta.gz",
    },
    "Cytb": {
        "NUC": "RAW/uniq/MIDORI2_UNIQ_NUC_GB272_Cytb_RAW.fasta.gz",
        "AA": "RAW_AA/total/MIDORI2_TOTAL_AA_GB272_Cytb_RAW_AA.fasta.gz",
    },
}

_SPECIES_TAXON = re.compile(r"(?:^|;)species_([^;]+)_([0-9]+)(?:;|$)")
_CLASS_TAXON = re.compile(r"(?:^|;)class_([^;]+)_([0-9]+)(?:;|$)")
_LOCUS_IDENTIFIER = re.compile(r"^(.+)\.([<>]?\d+)\.([<>]?\d+)$")
_AA_LOCUS_PREFIX = re.compile(r"^(.+\.[<>]?\d+\.[<>]?\d+)_")


def _normal_name(value: object) -> str:
    return " ".join(str(value).split()).casefold()


def _iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield identifier, taxonomy string, and sequence from gzipped FASTA."""

    identifier = ""
    taxonomy = ""
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if identifier:
                    yield identifier, taxonomy, "".join(chunks).upper()
                header = line[1:]
                parts = header.split(maxsplit=1)
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
    accession = re.sub(r"\.\d+$", "", accession_version)
    return {
        "accession_version": accession_version,
        "accession": accession,
        "start_text": start_text,
        "end_text": end_text,
        "is_partial": start_text.startswith("<") or end_text.startswith(">"),
    }


def _wanted_record(
    fields: Mapping[str, str],
    accepted_names: set[str],
    accepted_taxids: set[str],
) -> bool:
    return (
        fields["taxid"] in accepted_taxids
        or _normal_name(fields["organism"]) in accepted_names
    )


def _load_target_records(
    path: Path,
    accepted_names: set[str],
    accepted_taxids: set[str],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    records: list[dict[str, object]] = []
    counts = {"n_records": 0, "n_mammalia": 0, "n_target_candidates": 0}
    for identifier, taxonomy, sequence in _iter_fasta(path):
        counts["n_records"] += 1
        fields = _taxonomy_fields(taxonomy)
        if fields["class_name"] == "Mammalia" or fields["class_taxid"] == "40674":
            counts["n_mammalia"] += 1
        if not _wanted_record(fields, accepted_names, accepted_taxids):
            continue
        counts["n_target_candidates"] += 1
        records.append(
            {
                "identifier": identifier,
                "taxonomy": taxonomy,
                "sequence": sequence,
                **fields,
            }
        )
    return records, counts


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
            counts = dict((scan_counts or {}).get((gene, molecule), {}))
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
                    **counts,
                }
            )
    return pd.DataFrame(rows)


def _candidate_priority(
    record: Mapping[str, object],
    target: Mapping[str, object],
) -> int | None:
    taxid = str(record.get("taxid", ""))
    organism = _normal_name(record.get("organism", ""))
    resolved_taxid = str(target.get("Resolved_TaxID", ""))
    source_taxid = str(target.get("Source_TaxID", ""))
    original_name = _normal_name(target.get("Species_query", ""))
    resolved_name = _normal_name(target.get("Resolved_scientific_name", ""))
    if resolved_taxid and taxid == resolved_taxid:
        return 1
    if source_taxid and taxid == source_taxid:
        return 2
    if original_name and organism == original_name:
        return 3
    if resolved_name and organism == resolved_name:
        return 4
    return None


def _select_candidate(
    records: Iterable[Mapping[str, object]],
    target: Mapping[str, object],
) -> tuple[dict[str, object] | None, str]:
    ranked = [(_candidate_priority(record, target), dict(record)) for record in records]
    ranked = [(rank, record) for rank, record in ranked if rank is not None]
    if not ranked:
        return None, "not_found"
    best_rank = min(rank for rank, _ in ranked)
    best = [record for rank, record in ranked if rank == best_rank]
    unique = {str(record["identifier"]): record for record in best}
    if len(unique) != 1:
        return None, "ambiguous"
    return next(iter(unique.values())), "selected"


def _matching_protein(
    nucleotide: Mapping[str, object],
    proteins: Iterable[Mapping[str, object]],
) -> tuple[dict[str, object] | None, str]:
    prefix = f"{nucleotide['identifier']}_"
    paired = [
        dict(record)
        for record in proteins
        if str(record["identifier"]).startswith(prefix)
        and str(record.get("taxid", "")) == str(nucleotide.get("taxid", ""))
    ]
    unique = {str(record["identifier"]): record for record in paired}
    if len(unique) == 1:
        return next(iter(unique.values())), "same_accession_and_coordinates"
    if len(unique) > 1:
        return None, "ambiguous_same_locus_AA"
    return None, "no_same_locus_AA"


def _translation_parameters(cds: str, protein: str) -> tuple[int, bool, bool]:
    """Infer the retained CDS frame while requiring exact AA equality."""

    protein = protein.rstrip("*")
    for codon_start in (1, 2, 3):
        coding = cds[codon_start - 1 :]
        for apply_start_rule in (True, False):
            calculated = core.translate_table2(
                coding, apply_start_rule=apply_start_rule
            ).rstrip("*")
            if protein and calculated == protein:
                return codon_start, apply_start_rule, True
    return 1, True, False


_PROTEIN_LENGTH_BOUNDS = {
    "CO1": (450, 550),
    "CO3": (220, 300),
    "Cytb": (330, 420),
}


def _taxon_candidates(
    records: Iterable[Mapping[str, object]],
    target: Mapping[str, object],
) -> list[dict[str, object]]:
    ranked = [(_candidate_priority(record, target), dict(record)) for record in records]
    ranked = [(rank, record) for rank, record in ranked if rank is not None]
    if not ranked:
        return []
    best_rank = min(rank for rank, _ in ranked)
    return [record for rank, record in ranked if rank == best_rank]


def _prepare_candidate(
    species: str,
    gene: str,
    target: Mapping[str, object],
    nucleotide: Mapping[str, object],
    proteins: Iterable[Mapping[str, object]],
) -> dict[str, object]:
    protein, pair_status = _matching_protein(nucleotide, proteins)
    locus = _parse_locus_identifier(str(nucleotide["identifier"]))
    cds = str(nucleotide["sequence"])
    protein_sequence = str(protein["sequence"]) if protein else ""
    codon_start, apply_start_rule, inferred_match = _translation_parameters(
        cds, protein_sequence
    )
    coding = cds[codon_start - 1 :]
    remainder = len(coding) % 3
    terminal_fragment = coding[-remainder:] if remainder else ""
    frame_valid = remainder == 0 or terminal_fragment in {"T", "TA"}
    calculated = core.translate_table2(
        coding, apply_start_rule=apply_start_rule
    )
    internal_stop = "*" in calculated[:-1]
    ambiguous = sum(base not in "ACGT" for base in coding)
    ambiguous_fraction = ambiguous / len(coding) if coding else 1.0
    protein_length = len(protein_sequence.rstrip("*"))
    low, high = _PROTEIN_LENGTH_BOUNDS[gene]
    length_valid = low <= protein_length <= high
    basic_valid = bool(
        protein
        and pair_status == "same_accession_and_coordinates"
        and inferred_match
        and not bool(locus["is_partial"])
        and ambiguous_fraction <= 0.01
        and frame_valid
        and not internal_stop
        and length_valid
    )
    protein_accession = ""
    if protein:
        protein_accession = str(protein["identifier"])[
            len(str(nucleotide["identifier"])) + 1 :
        ]
    species_name_exact = _normal_name(nucleotide["organism"]) == _normal_name(
        target.get("Species_query", "")
    )
    resolved_taxid_exact = bool(str(target.get("Resolved_TaxID", ""))) and str(
        nucleotide["taxid"]
    ) == str(target.get("Resolved_TaxID", ""))
    match_status = "exact" if species_name_exact else (
        "taxid_confirmed" if resolved_taxid_exact else "manual_review"
    )
    is_refseq = str(locus["accession_version"]).startswith(("NC_", "NW_", "NZ_"))
    length_midpoint = (low + high) / 2
    sort_key = (
        0 if basic_valid else 1,
        0 if is_refseq else 1,
        0 if pair_status == "same_accession_and_coordinates" else 1,
        0 if inferred_match else 1,
        0 if not bool(locus["is_partial"]) else 1,
        ambiguous_fraction,
        abs(protein_length - length_midpoint),
        str(nucleotide["identifier"]),
    )
    sequence_row = {
        "Species": species,
        "Gene": gene,
        "display_gene": core.DISPLAY_GENE[gene],
        "accession": locus["accession"],
        "accession_version": locus["accession_version"],
        "TaxID": nucleotide["taxid"],
        "MIDORI_organism": nucleotide["organism"],
        "coding_strand": 1,
        "transl_table": 2,
        "cds_sequence": cds,
        "protein_sequence": protein_sequence,
        "cds_length": len(cds),
        "protein_length": protein_length,
        "gene_qualifier": "",
        "product_qualifier": "",
        "alias_match_method": "MIDORI2_gene_specific_database",
        "cds_coordinates": f"{locus['start_text']}..{locus['end_text']}",
        "codon_start": codon_start,
        "apply_start_rule": apply_start_rule,
        "is_partial": locus["is_partial"],
        "ambiguous_nucleotides": ambiguous,
        "has_ambiguous_nucleotides": bool(ambiguous),
        "record_source": "MIDORI2_UNIQ",
        "Source_TaxID": target.get("Source_TaxID", ""),
        "species_match_status": match_status,
        "species_match_method": (
            "normalized_exact_name" if species_name_exact else "resolved_NCBI_TaxID"
        ),
        "selection_reason": "",
        "midori_nuc_identifier": nucleotide["identifier"],
        "midori_aa_identifier": protein["identifier"] if protein else "",
        "protein_accession_version": protein_accession,
        "dna_aa_pair_status": pair_status,
        "same_feature_pair": pair_status == "same_accession_and_coordinates",
        "translation_match_inferred": inferred_match,
        "source_orientation": "coding_sequence_orientation",
        "midori_release": MIDORI_RELEASE,
    }
    manifest_row = {
        "Species": species,
        "Species_query": target.get("Species_query", ""),
        "Source_TaxID": target.get("Source_TaxID", ""),
        "Resolved_TaxID": target.get("Resolved_TaxID", ""),
        "Gene": gene,
        "selected": False,
        "selection_status": match_status,
        "midori_organism": nucleotide["organism"],
        "midori_taxid": nucleotide["taxid"],
        "nuc_identifier": nucleotide["identifier"],
        "aa_identifier": protein["identifier"] if protein else "",
        "accession_version": locus["accession_version"],
        "dna_aa_pair_status": pair_status,
        "is_refseq_accession": is_refseq,
        "is_partial": locus["is_partial"],
        "ambiguous_nucleotides": ambiguous,
        "ambiguous_fraction": ambiguous_fraction,
        "cds_length": len(cds),
        "protein_length": protein_length,
        "frame_valid": frame_valid,
        "internal_stop": internal_stop,
        "translation_match_inferred": inferred_match,
        "basic_valid": basic_valid,
        "release": MIDORI_RELEASE,
        "database_type": MIDORI_DATABASE_TYPE,
    }
    return {
        "accession_version": str(locus["accession_version"]),
        "basic_valid": basic_valid,
        "sort_key": sort_key,
        "sequence_row": sequence_row,
        "manifest_row": manifest_row,
        "nucleotide_record": dict(nucleotide),
        "match_status": match_status,
    }


def retrieve_midori_sequences(
    common_species: pd.DataFrame,
    cache_dir: str | Path,
    *,
    progress: Callable[[str], None] | None = print,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read MIDORI2 FASTAs and return sequence, audit, and match tables."""

    cache = Path(cache_dir)
    accepted_names = {
        _normal_name(value)
        for column in ("Species_query", "Resolved_scientific_name")
        if column in common_species
        for value in common_species[column]
        if str(value)
    }
    accepted_taxids = {
        str(value)
        for column in ("Source_TaxID", "Resolved_TaxID")
        if column in common_species
        for value in common_species[column]
        if str(value)
    }

    by_gene: dict[str, dict[str, list[dict[str, object]]]] = {}
    scan_counts: dict[tuple[str, str], dict[str, int]] = {}
    for gene in core.TARGET_GENES:
        by_gene[gene] = {}
        for molecule in ("NUC", "AA"):
            path = cache / Path(MIDORI_FILES[gene][molecule]).name
            if not path.exists():
                # Raises an error with the complete reproducible URL.
                midori_file_manifest(cache)
            if progress:
                progress(f"Scanning {path.name}")
            records, counts = _load_target_records(
                path, accepted_names, accepted_taxids
            )
            by_gene[gene][molecule] = records
            scan_counts[(gene, molecule)] = counts

    protein_index: dict[str, dict[str, list[dict[str, object]]]] = {}
    for gene in core.TARGET_GENES:
        index: dict[str, list[dict[str, object]]] = defaultdict(list)
        for record in by_gene[gene]["AA"]:
            match = _AA_LOCUS_PREFIX.match(str(record["identifier"]))
            if match:
                index[match.group(1)].append(record)
        protein_index[gene] = dict(index)

    sequences: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    match_rows: list[dict[str, object]] = []
    consistency_rows: list[dict[str, object]] = []
    for target in common_species.to_dict("records"):
        species = str(target["Species"])
        candidates_by_gene: dict[str, list[dict[str, object]]] = {}
        for gene in core.TARGET_GENES:
            nucleotide_candidates = _taxon_candidates(by_gene[gene]["NUC"], target)
            candidates_by_gene[gene] = sorted(
                [
                    _prepare_candidate(
                        species,
                        gene,
                        target,
                        nucleotide,
                        protein_index[gene].get(str(nucleotide["identifier"]), ()),
                    )
                    for nucleotide in nucleotide_candidates
                ],
                key=lambda candidate: candidate["sort_key"],
            )

        valid_accessions = {
            gene: {
                str(candidate["accession_version"])
                for candidate in candidates
                if candidate["basic_valid"] and candidate["accession_version"]
            }
            for gene, candidates in candidates_by_gene.items()
        }
        shared_accessions = set.intersection(
            *(valid_accessions[gene] for gene in core.TARGET_GENES)
        )
        if shared_accessions:
            def shared_score(accession: str) -> tuple[object, ...]:
                chosen = [
                    next(
                        candidate
                        for candidate in candidates_by_gene[gene]
                        if candidate["basic_valid"]
                        and candidate["accession_version"] == accession
                    )
                    for gene in core.TARGET_GENES
                ]
                return (
                    0 if accession.startswith(("NC_", "NW_", "NZ_")) else 1,
                    sum(float(candidate["sort_key"][5]) for candidate in chosen),
                    accession,
                )

            shared_accession = min(shared_accessions, key=shared_score)
            selected_by_gene = {
                gene: min(
                    (
                        candidate
                        for candidate in candidates_by_gene[gene]
                        if candidate["basic_valid"]
                        and candidate["accession_version"] == shared_accession
                    ),
                    key=lambda candidate: candidate["sort_key"],
                )
                for gene in core.TARGET_GENES
            }
            selection_mode = "shared_accession_all_three"
        else:
            selected_by_gene = {
                gene: (candidates[0] if candidates else None)
                for gene, candidates in candidates_by_gene.items()
            }
            selection_mode = "best_gene_specific_fallback"

        gene_statuses: dict[str, str] = {}
        for gene in core.TARGET_GENES:
            candidates = candidates_by_gene[gene]
            selected = selected_by_gene[gene]
            if selected is None:
                gene_statuses[gene] = "not_found"
                manifest_rows.append(
                    {
                        "Species": species,
                        "Species_query": target.get("Species_query", ""),
                        "Source_TaxID": target.get("Source_TaxID", ""),
                        "Resolved_TaxID": target.get("Resolved_TaxID", ""),
                        "Gene": gene,
                        "selected": False,
                        "candidate_rank": "",
                        "selection_status": "not_found",
                        "selection_reason": "no exact-name or resolved-TaxID MIDORI2 candidate",
                        "dna_aa_pair_status": "nucleotide_not_selected",
                        "release": MIDORI_RELEASE,
                        "database_type": MIDORI_DATABASE_TYPE,
                    }
                )
                continue
            if selection_mode == "shared_accession_all_three":
                selection_reason = (
                    "shared accession with basic-valid exact DNA-AA pairs for all three genes"
                )
            elif selected["basic_valid"]:
                selection_reason = (
                    "best basic-valid gene-specific exact DNA-AA pair; no valid shared accession"
                )
            else:
                selection_reason = (
                    "no basic-valid exact DNA-AA pair; best auditable candidate retained for QC"
                )
            selected["sequence_row"]["selection_reason"] = selection_reason
            sequences.append(selected["sequence_row"])
            gene_statuses[gene] = str(selected["match_status"])
            for rank, candidate in enumerate(candidates, start=1):
                manifest_row = dict(candidate["manifest_row"])
                is_selected = candidate is selected
                manifest_row.update(
                    candidate_rank=rank,
                    selected=is_selected,
                    selection_reason=selection_reason if is_selected else "not selected by deterministic score",
                    selection_mode=selection_mode,
                )
                manifest_rows.append(manifest_row)

        accessions = {
            gene: (str(candidate["accession_version"]) if candidate else "")
            for gene, candidate in selected_by_gene.items()
        }
        found_accessions = [value for value in accessions.values() if value]
        same_record = (
            len(found_accessions) == len(core.TARGET_GENES)
            and len(set(found_accessions)) == 1
        )
        reliable_statuses = {"exact", "taxid_confirmed"}
        statuses = list(gene_statuses.values())
        if any(status == "ambiguous" for status in statuses):
            overall_status = "ambiguous"
        elif not statuses or all(status == "not_found" for status in statuses):
            overall_status = "not_found"
        elif any(status == "manual_review" for status in statuses):
            overall_status = "manual_review"
        elif all(status in reliable_statuses or status == "not_found" for status in statuses):
            overall_status = (
                "exact" if all(status in {"exact", "not_found"} for status in statuses)
                else "taxid_confirmed"
            )
        else:
            overall_status = "manual_review"
        selected_records = [
            candidate["nucleotide_record"]
            for candidate in selected_by_gene.values()
            if candidate
        ]
        organisms = sorted({str(record["organism"]) for record in selected_records})
        taxids = sorted({str(record["taxid"]) for record in selected_records})
        selected_text = ";".join(
            f"{gene}:{accessions[gene]}" for gene in core.TARGET_GENES if accessions[gene]
        )
        match_rows.append(
            {
                "Species_original": species,
                "Species_query": target.get("Species_query", ""),
                "MIDORI_organism": ";".join(organisms),
                "TaxID": ";".join(taxids),
                "match_status": overall_status,
                "match_method": "exact_name_or_resolved_NCBI_TaxID",
                "selected_accession": selected_text,
                "notes": (
                    f"MIDORI2 {MIDORI_VERSION} UNIQ_NUC+TOTAL_AA; "
                    f"selection_mode={selection_mode}; same accession across all three={same_record}"
                ),
                "Source_TaxID": target.get("Source_TaxID", ""),
                "Resolved_TaxID": target.get("Resolved_TaxID", ""),
                "n_genes_found": len(found_accessions),
                "same_accession_all_three": same_record,
                "selection_mode": selection_mode,
            }
        )
        consistency_rows.append(
            {
                "Species": species,
                "CO1_accession_version": accessions["CO1"],
                "CO3_accession_version": accessions["CO3"],
                "Cytb_accession_version": accessions["Cytb"],
                "n_genes_found": len(found_accessions),
                "same_accession_all_three": same_record,
                "selection_mode": selection_mode,
                "all_selected_basic_valid": all(
                    candidate is not None and bool(candidate["basic_valid"])
                    for candidate in selected_by_gene.values()
                ),
            }
        )

    sequence_frame = pd.DataFrame(sequences)
    consistency = pd.DataFrame(consistency_rows)
    if not sequence_frame.empty:
        sequence_frame = sequence_frame.merge(
            consistency[["Species", "same_accession_all_three"]],
            on="Species",
            how="left",
            validate="many_to_one",
        )
    files = midori_file_manifest(cache, scan_counts)
    return (
        sequence_frame,
        pd.DataFrame(manifest_rows),
        pd.DataFrame(match_rows),
        consistency,
        files,
    )


def _write_downstream(
    spectra: pd.DataFrame,
    qc: pd.DataFrame,
    final_species: list[str],
    derived: Path,
    figures: Path,
    *,
    n_boot: int,
    seed: int,
) -> dict[str, object]:
    final_qc = qc.loc[qc["Species"].isin(final_species)].copy()
    composition, scores = core.observed_amino_acid_composition(final_qc)
    opportunities = core.codon_mutational_opportunities(final_qc)
    net = core.amino_acid_net_opportunities(opportunities)
    predicted = core.mutation_weighted_shifts(net, spectra, final_species)

    core._write_csv(composition, derived / "observed_amino_acid_composition.csv")
    core._write_csv(scores, derived / "observed_amino_acid_scores.csv")
    opportunities.to_csv(
        derived / "codon_mutational_opportunities.csv.gz",
        index=False,
        compression="gzip",
    )
    core._write_csv(net, derived / "amino_acid_net_opportunities.csv")
    core._write_csv(predicted, derived / "mutation_weighted_amino_acid_shifts.csv")

    score_long = scores.melt(
        ["Species", "Gene"], var_name="metric", value_name="value"
    )
    score_summary = core.cluster_bootstrap_mean(
        score_long, "value", ["Gene", "metric"], n_boot=n_boot, seed=seed
    )
    composition_summary = core.cluster_bootstrap_mean(
        composition,
        "aa_fraction",
        ["Gene", "amino_acid"],
        n_boot=n_boot,
        seed=seed + 1,
    )
    net_summary = core.cluster_bootstrap_mean(
        net,
        "normalized_net_opportunity",
        ["Gene", "mutation_H", "amino_acid"],
        n_boot=n_boot,
        seed=seed + 2,
    )
    gene_pairs = (
        ("COX3 - COX1", "CO3", "CO1"),
        ("CytB - COX3", "Cytb", "CO3"),
        ("CytB - COX1", "Cytb", "CO1"),
    )
    score_contrast_rows = []
    for metric in ("observed_A_rich_score", "observed_C_rich_score"):
        wide = scores.pivot(index="Species", columns="Gene", values=metric)
        for label, high, low in gene_pairs:
            difference = (wide[high] - wide[low]).rename("value").reset_index()
            difference["metric"] = metric
            difference["contrast"] = label
            score_contrast_rows.append(difference)
    score_contrast_data = pd.concat(score_contrast_rows, ignore_index=True)
    score_contrasts = core.cluster_bootstrap_mean(
        score_contrast_data,
        "value",
        ["metric", "contrast"],
        n_boot=n_boot,
        seed=seed + 3,
    )
    mutation_data = spectra.loc[
        spectra["Species"].isin(final_species)
        & spectra["Mut"].isin(core.MUTATION_TYPES),
        ["Species", "Gene", "Mut", "MutSpec"],
    ].drop_duplicates(["Species", "Gene", "Mut"])
    mutation_summary = core.cluster_bootstrap_mean(
        mutation_data,
        "MutSpec",
        ["Gene", "Mut"],
        n_boot=n_boot,
        seed=seed + 4,
    )
    predicted_summary = core.cluster_bootstrap_mean(
        predicted,
        "combined_predicted_shift",
        ["Gene", "amino_acid"],
        n_boot=n_boot,
        seed=seed + 5,
    )
    contrasts, associations = core.observed_predicted_contrasts(
        composition, predicted, n_boot=n_boot, seed=seed + 6
    )
    core._write_csv(score_summary, derived / "observed_amino_acid_score_summary.csv")
    core._write_csv(score_contrasts, derived / "observed_amino_acid_score_contrasts.csv")
    core._write_csv(
        composition_summary, derived / "observed_amino_acid_composition_summary.csv"
    )
    core._write_csv(net_summary, derived / "amino_acid_net_opportunity_summary.csv")
    core._write_csv(mutation_summary, derived / "mutation_pressure_gene_summary.csv")
    core._write_csv(predicted_summary, derived / "predicted_amino_acid_shift_summary.csv")
    core._write_csv(contrasts, derived / "observed_predicted_gene_contrasts.csv")
    core._write_csv(associations, derived / "observed_predicted_spearman.csv")
    core.plot_main_figure(
        score_summary,
        mutation_summary,
        predicted_summary,
        contrasts,
        associations,
        figures / "main_amino_acid_shift.png",
    )
    return {
        "composition": composition,
        "scores": scores,
        "opportunities": opportunities,
        "net_opportunities": net,
        "predicted": predicted,
        "score_summary": score_summary,
        "composition_summary": composition_summary,
        "net_summary": net_summary,
        "score_contrasts": score_contrasts,
        "mutation_summary": mutation_summary,
        "predicted_summary": predicted_summary,
        "contrasts": contrasts,
        "associations": associations,
    }


def run_midori_analysis(
    mutation_source: str | Path,
    stage_dir: str | Path,
    *,
    info_source: str | Path | None = None,
    class_filter: str | None = "Mammalia",
    n_boot: int = 2000,
    seed: int = 20260811,
    progress: Callable[[str], None] | None = print,
) -> dict[str, object]:
    """Execute stage 6 with MIDORI2 sequences and separate output folders."""

    if class_filter != "Mammalia":
        raise ValueError(
            "The MIDORI2 backend is currently restricted to Mammalia because "
            "sequence QC and codon consequences use vertebrate mitochondrial table 2"
        )

    stage = Path(stage_dir)
    cache = stage / "data" / "midori"
    derived = stage / "data" / "derived_midori"
    figures = stage / "figures" / "midori"
    for directory in (cache, derived, figures):
        directory.mkdir(parents=True, exist_ok=True)

    spectra, mutation_availability, common = core.load_mutation_cohort(
        mutation_source, class_filter=class_filter
    )
    mutation_gene_counts = (
        spectra[["Gene", "Species"]]
        .drop_duplicates()
        .groupby("Gene")
        .size()
        .reindex(core.TARGET_GENES)
        .rename("n_species")
        .reset_index()
    )
    core._write_csv(
        mutation_gene_counts, derived / "mutation_gene_species_counts.csv"
    )
    inferred_info = Path(mutation_source).parent / "info.csv"
    info_path = Path(info_source) if info_source is not None else inferred_info
    if info_path.exists():
        source_taxids = core.load_source_species_taxids(info_path)
        common = common.merge(
            source_taxids, on="Species_query", how="left", validate="one_to_one"
        )
        common["Source_TaxID"] = common["Source_TaxID"].fillna("").astype(str)
        taxid_map = core.NCBITaxonomyClient(cache).resolve_taxids(
            common.loc[common["Source_TaxID"].ne(""), "Source_TaxID"].tolist()
        )
        common = common.merge(
            taxid_map, on="Source_TaxID", how="left", validate="many_to_one"
        )
    else:
        common["Source_TaxID"] = ""
        common["source_taxid_status"] = "source_info_unavailable"
        common["Resolved_TaxID"] = ""
        common["Resolved_scientific_name"] = ""
        common["TaxID_was_merged"] = False
    for column in ("Resolved_TaxID", "Resolved_scientific_name"):
        common[column] = common[column].fillna("").astype(str)
    common["TaxID_was_merged"] = common["TaxID_was_merged"].fillna(False).astype(bool)
    species_taxids = common[
        [
            "Species",
            "Species_query",
            "Source_TaxID",
            "source_taxid_status",
            "Resolved_TaxID",
            "Resolved_scientific_name",
            "TaxID_was_merged",
        ]
    ].copy()
    core._write_csv(species_taxids, derived / "mutation_species_taxids.csv")

    sequences, manifest, matches, consistency, file_manifest = (
        retrieve_midori_sequences(common, cache, progress=progress)
    )
    qc = core.sequence_qc(sequences)
    availability, attrition, exclusions = core.build_availability_and_attrition(
        mutation_availability,
        common,
        matches,
        qc,
        source_label="MIDORI2",
    )
    availability = availability.merge(
        consistency, on="Species", how="left", validate="one_to_one"
    )
    exclusion_order = (
        "not found in MIDORI2",
        "ambiguous taxonomic match",
        "gene missing",
        "partial/problematic CDS",
        "translation mismatch",
        "other QC failure",
    )
    exclusion_summary = (
        exclusions.groupby("exclusion_reason")
        .size()
        .reindex(exclusion_order, fill_value=0)
        .rename("n_species")
        .reset_index()
    )
    final_species = sorted(
        availability.loc[availability["included_final"], "Species"]
    )

    core._write_csv(file_manifest, derived / "midori_file_manifest.csv")
    core._write_csv(manifest, derived / "midori_sequence_manifest.csv")
    core._write_csv(matches, derived / "species_midori_matching.csv")
    core._write_csv(consistency, derived / "midori_record_consistency.csv")
    core._write_csv(qc, derived / "sequence_qc.csv")
    core._write_csv(availability, derived / "species_gene_availability.csv")
    core._write_csv(attrition, derived / "species_intersection_summary.csv")
    core._write_csv(
        pd.DataFrame({"Species": final_species}),
        derived / "final_common_species.csv",
    )
    core._write_csv(exclusions, derived / "species_exclusion_reasons.csv")
    core._write_csv(exclusion_summary, derived / "exclusion_reason_summary.csv")
    core._write_csv(
        core.amino_acid_codon_richness(),
        derived / "amino_acid_codon_richness.csv",
    )
    canonical_first = [
        "Species",
        "Gene",
        "display_gene",
        "accession",
        "accession_version",
        "TaxID",
        "MIDORI_organism",
        "coding_strand",
        "transl_table",
        "cds_sequence",
        "protein_sequence",
        "cds_length",
        "protein_length",
        "sequence_qc_status",
    ]
    canonical = qc[
        [
            *canonical_first,
            *[column for column in qc.columns if column not in canonical_first],
        ]
    ]
    canonical.to_csv(
        derived / "midori_gene_sequences.csv.gz", index=False, compression="gzip"
    )
    outputs: dict[str, object] = {
        "spectra": spectra,
        "mutation_gene_counts": mutation_gene_counts,
        "species_taxids": species_taxids,
        "common_species": common,
        "sequences": sequences,
        "manifest": manifest,
        "matches": matches,
        "record_consistency": consistency,
        "file_manifest": file_manifest,
        "qc": qc,
        "availability": availability,
        "attrition": attrition,
        "exclusions": exclusions,
        "exclusion_summary": exclusion_summary,
        "final_species": final_species,
        "retrieval_backend": f"MIDORI2 {MIDORI_VERSION} UNIQ_NUC+TOTAL_AA RAW",
        "source_limitation": (
            "UNIQ_NUC collapses identical nucleotide haplotypes; a shared accession "
            "is preferred when represented, and same_accession_all_three is audited"
        ),
    }
    if final_species:
        outputs.update(
            _write_downstream(
                spectra,
                qc,
                final_species,
                derived,
                figures,
                n_boot=n_boot,
                seed=seed,
            )
        )
    return outputs


__all__ = [
    "MIDORI_BASE_URL",
    "MIDORI_DATABASE_TYPE",
    "MIDORI_FILES",
    "MIDORI_RELEASE",
    "MIDORI_VERSION",
    "midori_file_manifest",
    "retrieve_midori_sequences",
    "run_midori_analysis",
]
