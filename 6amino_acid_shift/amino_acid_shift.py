"""Shared calculations for the MIDORI2 mitochondrial amino-acid analysis."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from io import BytesIO
from http.client import IncompleteRead
import json
import os
from pathlib import Path
import time
from typing import Callable, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:  # Optional by design: the task explicitly forbids automatic installs.
    from Bio import Entrez

    BIOPYTHON_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised in minimal environments.
    Entrez = None
    BIOPYTHON_AVAILABLE = False


TARGET_GENES = ("CO1", "CO3", "Cytb")
DISPLAY_GENE = {"CO1": "COX1", "CO3": "COX3", "Cytb": "CytB"}
GENE_ORDER = {gene: index for index, gene in enumerate(TARGET_GENES)}
MUTATION_TYPES = ("A>G", "C>T")
MUTATION_DISPLAY = {"A>G": "A_H>G_H", "C>T": "C_H>T_H"}
AA_ORDER = tuple("ACDEFGHIKLMNPQRSTVWY")
# Vertebrate mitochondrial code (NCBI translation table 2). Stops are '*'.
_TABLE2_ROWS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "*", "AGG": "*",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_TABLE2_STARTS = frozenset(("ATT", "ATC", "ATA", "ATG", "GTG"))


def safe_species_query(species_identifier: object) -> str:
    """Convert repository underscores to spaces without dropping taxon tokens."""

    return " ".join(str(species_identifier).replace("_", " ").split())


def load_mutation_cohort(
    source: str | Path,
    *,
    class_filter: str | None = "Mammalia",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load spectra and return target rows, availability, and exact intersection."""

    spectra = pd.read_csv(source)
    required = {"Species", "Gene", "Mut", "MutSpec"}
    missing = required - set(spectra.columns)
    if missing:
        raise ValueError(f"Mutation-spectrum input lacks columns: {sorted(missing)}")
    spectra = spectra[spectra["Gene"].isin(TARGET_GENES)].copy()
    if class_filter is not None:
        if "Class" not in spectra.columns:
            raise ValueError("class_filter was supplied but the source has no Class column")
        spectra = spectra[spectra["Class"].eq(class_filter)].copy()
    spectra["Species"] = spectra["Species"].astype(str)

    presence = (
        spectra[["Species", "Gene"]]
        .drop_duplicates()
        .assign(present=True)
        .pivot(index="Species", columns="Gene", values="present")
        .notna()
    )
    for gene in TARGET_GENES:
        if gene not in presence:
            presence[gene] = False
    availability = presence[list(TARGET_GENES)].rename(
        columns={gene: f"has_mut_{gene}" for gene in TARGET_GENES}
    ).reset_index()
    common = availability.loc[
        availability[[f"has_mut_{gene}" for gene in TARGET_GENES]].all(axis=1),
        ["Species"],
    ].copy()
    common["Species_query"] = common["Species"].map(safe_species_query)
    common = common.sort_values("Species", kind="stable").reset_index(drop=True)
    return spectra, availability, common


def load_source_species_taxids(source: str | Path) -> pd.DataFrame:
    """Extract an unambiguous terminal species TaxID from the source taxonomy."""

    info = pd.read_csv(source, usecols=["species", "taxa"], dtype=str).drop_duplicates()
    info["Species_query"] = info["species"].map(safe_species_query)
    info["Source_TaxID"] = info["taxa"].str.extract(r"_(\d+)$", expand=False)
    rows = []
    for species_query, group in info.groupby("Species_query", sort=False):
        taxids = sorted(set(group["Source_TaxID"].dropna()))
        rows.append(
            {
                "Species_query": species_query,
                "Source_TaxID": taxids[0] if len(taxids) == 1 else "",
                "source_taxid_status": "unique" if len(taxids) == 1 else (
                    "missing" if not taxids else "ambiguous"
                ),
            }
        )
    return pd.DataFrame(rows)


@dataclass
class RateLimiter:
    interval_seconds: float
    last_request: float = 0.0

    def wait(self) -> None:
        delay = self.interval_seconds - (time.monotonic() - self.last_request)
        if delay > 0:
            time.sleep(delay)
        self.last_request = time.monotonic()


@dataclass(frozen=True)
class EntrezSettings:
    email: str | None
    api_key: str | None
    tool: str = "TBSSSignature"
    max_attempts: int = 4

    @classmethod
    def from_environment(cls) -> "EntrezSettings":
        return cls(
            email=os.getenv("NCBI_EMAIL") or os.getenv("ENTREZ_EMAIL"),
            api_key=os.getenv("NCBI_API_KEY") or os.getenv("ENTREZ_API_KEY"),
        )


class NCBITaxonomyClient:
    """Resolve current NCBI TaxIDs, with a small local mapping cache."""

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, cache_dir: str | Path, settings: EntrezSettings | None = None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.settings = settings or EntrezSettings.from_environment()
        self.backend = "Biopython" if BIOPYTHON_AVAILABLE else "stdlib_taxonomy_xml"
        # NCBI allows 10 requests/s with a key and 3 requests/s without one.
        self.limiter = RateLimiter(0.11 if self.settings.api_key else 0.35)
        if BIOPYTHON_AVAILABLE:
            Entrez.email = self.settings.email
            Entrez.api_key = self.settings.api_key
            Entrez.tool = self.settings.tool

    def _retry(self, operation: Callable[[], str]) -> str:
        last_error: Exception | None = None
        for attempt in range(self.settings.max_attempts):
            try:
                self.limiter.wait()
                return operation()
            except (HTTPError, URLError, OSError, IncompleteRead, TimeoutError) as error:
                last_error = error
                code = getattr(error, "code", None)
                if code is not None and code < 500 and code != 429:
                    raise
                if attempt + 1 < self.settings.max_attempts:
                    time.sleep(2**attempt)
        raise RuntimeError("NCBI request failed after retries") from last_error

    def _stdlib_get(self, endpoint: str, parameters: Mapping[str, object]) -> str:
        params = {**parameters, "tool": self.settings.tool}
        if self.settings.email:
            params["email"] = self.settings.email
        if self.settings.api_key:
            params["api_key"] = self.settings.api_key
        url = f"{self.base_url}/{endpoint}?{urlencode(params)}"

        def request() -> str:
            with urlopen(Request(url, headers={"User-Agent": self.settings.tool}), timeout=60) as response:
                return response.read().decode("utf-8")

        return self._retry(request)

    def resolve_taxids(
        self,
        source_taxids: Sequence[str],
        *,
        force_refresh: bool = False,
    ) -> pd.DataFrame:
        """Resolve current/merged NCBI TaxIDs and cache the mapping."""

        cache_path = self.cache_dir / "taxonomy_taxid_map.json"
        cached: dict[str, dict[str, object]] = {}
        if cache_path.exists() and not force_refresh:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
        requested = list(dict.fromkeys(str(value) for value in source_taxids if str(value)))
        missing = [value for value in requested if value not in cached]
        if missing:
            if BIOPYTHON_AVAILABLE:
                def request() -> bytes:
                    handle = Entrez.efetch(db="taxonomy", id=",".join(missing), retmode="xml")
                    try:
                        return handle.read()
                    finally:
                        handle.close()

                payload = self._retry(request)
                records = Entrez.read(BytesIO(payload))
                parsed = [
                    {
                        "current_taxid": str(record.get("TaxId", "")),
                        "current_scientific_name": str(record.get("ScientificName", "")),
                        "aka_taxids": [str(value) for value in record.get("AkaTaxIds", ())],
                    }
                    for record in records
                ]
            else:
                payload = self._stdlib_get(
                    "efetch.fcgi",
                    {"db": "taxonomy", "id": ",".join(missing), "retmode": "xml"},
                )
                root = ET.fromstring(payload)
                parsed = []
                for record in root.findall("./Taxon"):
                    parsed.append(
                        {
                            "current_taxid": record.findtext("TaxId", default=""),
                            "current_scientific_name": record.findtext("ScientificName", default=""),
                            "aka_taxids": [node.text or "" for node in record.findall("./AkaTaxIds/TaxId")],
                        }
                    )
            for source_taxid in missing:
                match = next(
                    (
                        record for record in parsed
                        if source_taxid == record["current_taxid"] or source_taxid in record["aka_taxids"]
                    ),
                    None,
                )
                cached[source_taxid] = match or {
                    "current_taxid": "",
                    "current_scientific_name": "",
                    "aka_taxids": [],
                }
            cache_path.write_text(json.dumps(cached, indent=2, sort_keys=True), encoding="utf-8")
        return pd.DataFrame(
            [
                {
                    "Source_TaxID": source_taxid,
                    "Resolved_TaxID": cached.get(source_taxid, {}).get("current_taxid", ""),
                    "Resolved_scientific_name": cached.get(source_taxid, {}).get("current_scientific_name", ""),
                    "TaxID_was_merged": source_taxid in cached.get(source_taxid, {}).get("aka_taxids", ()),
                }
                for source_taxid in requested
            ]
        )

def translate_table2(coding_sequence: str, *, apply_start_rule: bool = True) -> str:
    codons = [coding_sequence[i:i + 3] for i in range(0, len(coding_sequence) - 2, 3)]
    protein = [_TABLE2_ROWS.get(codon, "X") for codon in codons]
    if apply_start_rule and codons and codons[0] in _TABLE2_STARTS:
        protein[0] = "M"
    return "".join(protein)


def sequence_qc(sequence_data: pd.DataFrame, *, maximum_ambiguous_fraction: float = 0.01) -> pd.DataFrame:
    """Apply explicit per-feature QC, including terminal incomplete stops."""

    if sequence_data.empty:
        columns = list(sequence_data.columns) + [
            "coding_cds_sequence", "computed_translation", "translation_match",
            "terminal_incomplete_stop", "sequence_qc_status", "sequence_qc_reason",
        ]
        return pd.DataFrame(columns=list(dict.fromkeys(columns)))
    rows = []
    for source in sequence_data.to_dict("records"):
        row = dict(source)
        reasons: list[str] = []
        cds = str(row.get("cds_sequence", "")).upper()
        protein = str(row.get("protein_sequence", "")).replace(" ", "").upper().rstrip("*")
        codon_start = int(row.get("codon_start", 1) or 1)
        coding = cds[codon_start - 1:] if codon_start in (1, 2, 3) else cds
        table = row.get("transl_table")
        if not cds:
            reasons.append("missing_or_empty_CDS")
        if row.get("coding_strand") not in (-1, 1):
            reasons.append("coding_strand_undefined")
        if pd.isna(table) or int(table) != 2:
            reasons.append("unexpected_or_missing_transl_table")
        if codon_start not in (1, 2, 3):
            reasons.append("invalid_codon_start")
        if row.get("is_partial"):
            reasons.append("partial_CDS")
        ambiguous = sum(base not in "ACGT" for base in coding)
        if coding and ambiguous / len(coding) > maximum_ambiguous_fraction:
            reasons.append("excess_ambiguous_nucleotides")
        remainder = len(coding) % 3
        terminal_fragment = coding[-remainder:] if remainder else ""
        terminal_incomplete = terminal_fragment in {"T", "TA"}
        if remainder and not terminal_incomplete:
            reasons.append("frameshift_or_unrecognized_terminal_fragment")
        # Mitochondrial CDS features use the alternative start-codon rule.
        # MIDORI2 can also contain partial gene records, for which the
        # first retained codon is not necessarily the biological start codon.
        # A source adapter may therefore disable only that first-codon rule;
        # all remaining codons still use vertebrate mitochondrial table 2.
        apply_start_rule = bool(row.get("apply_start_rule", True))
        calculated = translate_table2(coding, apply_start_rule=apply_start_rule)
        if "*" in calculated[:-1]:
            reasons.append("internal_stop_codon")
        calculated_protein = calculated[:-1] if calculated.endswith("*") else calculated
        if not protein:
            reasons.append("missing_protein_translation")
        translation_match = bool(protein) and calculated_protein == protein
        if protein and not translation_match:
            reasons.append("translation_mismatch")
        row.update(
            coding_cds_sequence=coding,
            computed_translation=calculated_protein,
            translation_match=translation_match,
            terminal_incomplete_stop=terminal_incomplete,
            ambiguous_nucleotides=ambiguous,
        )
        rows.append((row, reasons))

    result = pd.DataFrame([row for row, _ in rows])
    # Length plausibility is relative to the retrieved distribution for the same gene.
    fallback_bounds = {"CO1": (450, 550), "CO3": (220, 300), "Cytb": (330, 420)}
    gene_bounds: dict[str, tuple[float, float]] = {}
    for gene, group in result.groupby("Gene"):
        lengths = group.loc[group["protein_length"].gt(0), "protein_length"].astype(float)
        if len(lengths) >= 3:
            median = float(lengths.median())
            mad = float((lengths - median).abs().median())
            tolerance = max(0.15 * median, 5.0 * mad, 5.0)
            gene_bounds[gene] = (median - tolerance, median + tolerance)
        else:
            gene_bounds[gene] = fallback_bounds[gene]
    for index, (row, reasons) in enumerate(rows):
        low, high = gene_bounds[str(row["Gene"])]
        length = float(row.get("protein_length", 0))
        if not (low <= length <= high):
            reasons.append(f"protein_length_outlier[{low:.1f},{high:.1f}]")
        result.loc[index, "sequence_qc_status"] = "pass" if not reasons else "fail"
        result.loc[index, "sequence_qc_reason"] = "pass" if not reasons else ";".join(reasons)
    return result


def build_availability_and_attrition(
    mutation_availability: pd.DataFrame,
    common_species: pd.DataFrame,
    matches: pd.DataFrame,
    qc: pd.DataFrame,
    *,
    source_label: str = "MIDORI2",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = common_species[["Species"]].merge(mutation_availability, on="Species", how="left")
    reliable = set(
        matches.loc[matches["match_status"].isin(("exact", "taxid_confirmed")), "Species_original"]
    ) if not matches.empty else set()
    for gene in TARGET_GENES:
        seq_species = set(qc.loc[qc["Gene"].eq(gene), "Species"]) if not qc.empty else set()
        valid_species = set(
            qc.loc[qc["Gene"].eq(gene) & qc["sequence_qc_status"].eq("pass"), "Species"]
        ) if not qc.empty else set()
        base[f"has_seq_{gene}"] = base["Species"].isin(seq_species) & base["Species"].isin(reliable)
        base[f"valid_{gene}"] = base["Species"].isin(valid_species) & base["Species"].isin(reliable)
    base["included_final"] = base[[f"valid_{gene}" for gene in TARGET_GENES]].all(axis=1)
    base = base[
        [
            "Species",
            *[f"has_mut_{gene}" for gene in TARGET_GENES],
            *[f"has_seq_{gene}" for gene in TARGET_GENES],
            *[f"valid_{gene}" for gene in TARGET_GENES],
            "included_final",
        ]
    ]

    n_mut = len(base)
    n_match = int(base["Species"].isin(reliable).sum())
    n_all_seq = int(base[[f"has_seq_{gene}" for gene in TARGET_GENES]].all(axis=1).sum())
    n_qc = int(base["included_final"].sum())
    counts = [n_mut, n_match, n_all_seq, n_qc, n_qc]
    stages = [
        ("mutation-spectrum cohort", "complete spectra for CO1, CO3 and Cytb"),
        (f"{source_label} matched", "not found or unsafe taxonomic match"),
        ("all 3 genes found", "one or more target CDS absent"),
        ("sequence QC", "one or more target CDS failed sequence QC"),
        ("final analysis cohort", "complete matched analysis input"),
    ]
    summary = []
    previous = counts[0]
    for index, ((stage, reason), count) in enumerate(zip(stages, counts)):
        summary.append(
            {"stage": stage, "n_species": count, "n_excluded": 0 if index == 0 else previous - count, "reason": reason}
        )
        previous = count

    match_lookup = matches.set_index("Species_original") if not matches.empty else pd.DataFrame()
    exclusions = []
    for item in base.loc[~base["included_final"]].itertuples(index=False):
        species = item.Species
        if species in getattr(match_lookup, "index", ()):
            match_row = match_lookup.loc[species]
            if isinstance(match_row, pd.DataFrame):
                match_row = match_row.iloc[0]
            match_status = match_row.get("match_status", "not_found")
            match_notes = str(match_row.get("notes", ""))
            selected_accession = str(match_row.get("selected_accession", ""))
        else:
            match_status, match_notes, selected_accession = "not_found", "No match row", ""
        missing_genes = [gene for gene in TARGET_GENES if not getattr(item, f"has_seq_{gene}")]
        failed = qc.loc[qc["Species"].eq(species) & qc["sequence_qc_status"].eq("fail")]
        qc_text = ";".join(failed.get("sequence_qc_reason", pd.Series(dtype=str)).astype(str))
        if match_status == "not_found":
            reason = f"not found in {source_label}"
        elif match_status in {"manual_review", "ambiguous"}:
            reason = "ambiguous taxonomic match"
        elif missing_genes:
            reason = "gene missing"
        else:
            if "partial_CDS" in qc_text or "frameshift" in qc_text or "protein_length_outlier" in qc_text:
                reason = "partial/problematic CDS"
            elif "translation_mismatch" in qc_text:
                reason = "translation mismatch"
            else:
                reason = "other QC failure"
        exclusions.append(
            {
                "Species": species,
                "exclusion_reason": reason,
                "match_status": match_status,
                "missing_genes": ";".join(missing_genes),
                "sequence_qc_details": qc_text,
                "selected_accession": selected_accession,
                "notes": match_notes,
            }
        )
    return base, pd.DataFrame(summary), pd.DataFrame(exclusions)


def amino_acid_codon_richness() -> pd.DataFrame:
    """Mean base fraction among synonymous table-2 codons for each amino acid."""

    rows = []
    for aa in AA_ORDER:
        codons = sorted(codon for codon, encoded in _TABLE2_ROWS.items() if encoded == aa)
        rows.append(
            {
                "amino_acid": aa,
                "n_synonymous_codons_table2": len(codons),
                "synonymous_codons_table2": ";".join(codons),
                "A_richness_12": np.mean([codon.count("A") / 3 for codon in codons]),
                "C_richness_12": np.mean([codon.count("C") / 3 for codon in codons]),
            }
        )
    return pd.DataFrame(rows)


def observed_amino_acid_composition(qc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    valid = qc.loc[qc["sequence_qc_status"].eq("pass")].copy()
    richness = amino_acid_codon_richness().set_index("amino_acid")
    composition, scores = [], []
    for row in valid.itertuples(index=False):
        protein = str(row.protein_sequence).rstrip("*")
        counts = Counter(protein)
        denominator = sum(counts.get(aa, 0) for aa in AA_ORDER)
        fractions = {aa: counts.get(aa, 0) / denominator for aa in AA_ORDER}
        for aa in AA_ORDER:
            composition.append(
                {
                    "Species": row.Species, "Gene": row.Gene,
                    "amino_acid": aa, "aa_count": counts.get(aa, 0),
                    "protein_length": denominator, "aa_fraction": fractions[aa],
                }
            )
        scores.append(
            {
                "Species": row.Species,
                "Gene": row.Gene,
                "observed_A_rich_score": sum(fractions[aa] * richness.loc[aa, "A_richness_12"] for aa in AA_ORDER),
                "observed_C_rich_score": sum(fractions[aa] * richness.loc[aa, "C_richness_12"] for aa in AA_ORDER),
            }
        )
    return pd.DataFrame(composition), pd.DataFrame(scores)


def codon_mutational_opportunities(qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    valid = qc.loc[qc["sequence_qc_status"].eq("pass")]
    for item in valid.itertuples(index=False):
        coding = str(item.coding_cds_sequence)
        codons = [coding[index:index + 3] for index in range(0, len(coding) - 2, 3)]
        for codon_index, ref_codon in enumerate(codons, start=1):
            if set(ref_codon) - set("ACGT"):
                continue
            ref_aa = _TABLE2_ROWS[ref_codon]
            for codon_position, base in enumerate(ref_codon, start=1):
                alternate = {"A": "G", "C": "T"}.get(base)
                if alternate is None:
                    continue
                alt_codon = ref_codon[:codon_position - 1] + alternate + ref_codon[codon_position:]
                alt_aa = _TABLE2_ROWS[alt_codon]
                if ref_aa == "*" and alt_aa != "*":
                    consequence = "stop_loss"
                elif ref_aa != "*" and alt_aa == "*":
                    consequence = "nonsense"
                elif ref_aa == alt_aa:
                    consequence = "synonymous"
                else:
                    consequence = "missense"
                rows.append(
                    {
                        "Species": item.Species, "Gene": item.Gene,
                        "codon_index": codon_index, "codon_position": codon_position,
                        "ref_codon": ref_codon, "alt_codon": alt_codon,
                        "ref_aa": ref_aa, "alt_aa": alt_aa,
                        "mutation_H": MUTATION_DISPLAY[f"{base}>{alternate}"],
                        "mutation": f"{base}>{alternate}", "consequence": consequence,
                    }
                )
    return pd.DataFrame(rows)


def amino_acid_net_opportunities(opportunities: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if opportunities.empty:
        return pd.DataFrame(columns=[
            "Species", "Gene", "mutation_H", "mutation", "amino_acid",
            "outgoing_opportunities", "incoming_opportunities", "net_opportunity",
            "total_valid_opportunities", "normalized_net_opportunity",
        ])
    for keys, group in opportunities.groupby(["Species", "Gene", "mutation_H", "mutation"], sort=False):
        total = len(group)
        changed = group.loc[group["ref_aa"].ne(group["alt_aa"])]
        for aa in AA_ORDER:
            outgoing = int(changed["ref_aa"].eq(aa).sum())
            incoming = int(changed["alt_aa"].eq(aa).sum())
            rows.append(
                {
                    "Species": keys[0], "Gene": keys[1], "mutation_H": keys[2],
                    "mutation": keys[3], "amino_acid": aa,
                    "outgoing_opportunities": outgoing,
                    "incoming_opportunities": incoming,
                    "net_opportunity": incoming - outgoing,
                    "total_valid_opportunities": total,
                    "normalized_net_opportunity": (incoming - outgoing) / total if total else np.nan,
                }
            )
    return pd.DataFrame(rows)


def mutation_weighted_shifts(net: pd.DataFrame, spectra: pd.DataFrame, final_species: Sequence[str]) -> pd.DataFrame:
    weights = spectra.loc[
        spectra["Species"].isin(final_species) & spectra["Mut"].isin(MUTATION_TYPES),
        ["Species", "Gene", "Mut", "MutSpec"],
    ].copy()
    duplicates = weights.duplicated(["Species", "Gene", "Mut"], keep=False)
    if duplicates.any():
        varying = weights.loc[duplicates].groupby(["Species", "Gene", "Mut"])["MutSpec"].nunique()
        if varying.gt(1).any():
            raise ValueError("Mutation input has conflicting duplicate MutSpec values")
        weights = weights.drop_duplicates(["Species", "Gene", "Mut"])
    joined = net.merge(
        weights.rename(columns={"Mut": "mutation"}),
        on=["Species", "Gene", "mutation"],
        how="left",
        validate="many_to_one",
    )
    if joined["MutSpec"].isna().any():
        raise ValueError("At least one final Species × Gene lacks A>G or C>T MutSpec")
    joined["weighted_contribution"] = joined["MutSpec"] * joined["normalized_net_opportunity"]
    wide = joined.pivot_table(
        index=["Species", "Gene", "amino_acid"], columns="mutation",
        values="weighted_contribution", aggfunc="first",
    ).reset_index()
    for mutation in MUTATION_TYPES:
        if mutation not in wide:
            wide[mutation] = 0.0
    wide = wide.rename(columns={"A>G": "AG_contribution", "C>T": "CT_contribution"})
    wide["combined_predicted_shift"] = wide["AG_contribution"] + wide["CT_contribution"]
    return wide


def cluster_bootstrap_mean(
    data: pd.DataFrame,
    value_col: str,
    group_cols: Sequence[str],
    *,
    n_boot: int = 2000,
    seed: int = 20260811,
) -> pd.DataFrame:
    species = np.array(sorted(pd.unique(data["Species"])), dtype=object)
    collapsed = data.groupby(["Species", *group_cols], sort=False)[value_col].mean().reset_index()
    observed = collapsed.groupby(list(group_cols), sort=False)[value_col].mean().rename("mean")
    if len(species) == 0:
        return observed.reset_index().assign(ci_low=np.nan, ci_high=np.nan, n_species=0)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(species), size=(n_boot, len(species)))
    replicate_rows = []
    for group_key, group in collapsed.groupby(list(group_cols), sort=False):
        values = group.set_index("Species")[value_col].reindex(species).to_numpy(float)
        boot_values = np.nanmean(values[draws], axis=1)
        key_tuple = group_key if isinstance(group_key, tuple) else (group_key,)
        replicate_rows.append((*key_tuple, np.quantile(boot_values, 0.025), np.quantile(boot_values, 0.975)))
    boot = pd.DataFrame(replicate_rows, columns=[*group_cols, "ci_low", "ci_high"]).set_index(list(group_cols))
    result = observed.to_frame()
    result = result.join(boot)
    result["n_species"] = len(species)
    return result.reset_index()


def observed_predicted_contrasts(
    composition: pd.DataFrame,
    predicted: pd.DataFrame,
    *,
    n_boot: int = 2000,
    seed: int = 20260811,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_wide = composition.pivot(index=["Species", "amino_acid"], columns="Gene", values="aa_fraction")
    predicted_wide = predicted.pivot(index=["Species", "amino_acid"], columns="Gene", values="combined_predicted_shift")
    pairs = (("COX3 - COX1", "CO3", "CO1"), ("CytB - COX3", "Cytb", "CO3"), ("CytB - COX1", "Cytb", "CO1"))
    contrast_rows = []
    for label, high, low in pairs:
        joined = pd.DataFrame(
            {
                "observed": observed_wide[high] - observed_wide[low],
                "predicted": predicted_wide[high] - predicted_wide[low],
            }
        ).dropna().reset_index()
        means = joined.groupby("amino_acid")[["observed", "predicted"]].mean().reset_index()
        means.insert(0, "contrast", label)
        contrast_rows.append(means)
    contrasts = pd.concat(contrast_rows, ignore_index=True)

    species = np.array(sorted(pd.unique(composition["Species"])), dtype=object)
    amino_acids = np.array(AA_ORDER, dtype=object)
    genes = np.array(TARGET_GENES, dtype=object)
    observed_cube = (
        composition.set_index(["Species", "amino_acid", "Gene"])["aa_fraction"]
        .reindex(pd.MultiIndex.from_product([species, amino_acids, genes]))
        .to_numpy(float).reshape(len(species), len(amino_acids), len(genes))
    )
    predicted_cube = (
        predicted.set_index(["Species", "amino_acid", "Gene"])["combined_predicted_shift"]
        .reindex(pd.MultiIndex.from_product([species, amino_acids, genes]))
        .to_numpy(float).reshape(len(species), len(amino_acids), len(genes))
    )
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(species), size=(n_boot, len(species)))
    observed_boot = np.nanmean(observed_cube[draws], axis=1)
    predicted_boot = np.nanmean(predicted_cube[draws], axis=1)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    association_rows = []
    for label, high, low in pairs:
        table = contrasts.loc[contrasts["contrast"].eq(label)]
        rho = float(spearmanr(table["predicted"], table["observed"]).statistic)
        high_index, low_index = gene_index[high], gene_index[low]
        boot_rho = [
            float(spearmanr(
                predicted_boot[replicate, :, high_index] - predicted_boot[replicate, :, low_index],
                observed_boot[replicate, :, high_index] - observed_boot[replicate, :, low_index],
            ).statistic)
            for replicate in range(n_boot)
        ]
        finite = np.array([value for value in boot_rho if np.isfinite(value)])
        association_rows.append(
            {
                "contrast": label, "spearman_rho": rho,
                "ci_low": np.quantile(finite, 0.025) if len(finite) else np.nan,
                "ci_high": np.quantile(finite, 0.975) if len(finite) else np.nan,
                "n_amino_acids": len(table), "n_species": len(species),
                "bootstrap_unit": "species cluster",
            }
        )
    return contrasts, pd.DataFrame(association_rows)


def _errorbar(ax: object, summary: pd.DataFrame, value_name: str, label: str, color: str) -> None:
    ordered = summary.set_index("Gene").reindex(TARGET_GENES)
    x = np.arange(len(TARGET_GENES))
    y = ordered["mean"].to_numpy(float)
    err = np.vstack((y - ordered["ci_low"].to_numpy(float), ordered["ci_high"].to_numpy(float) - y))
    ax.errorbar(x, y, yerr=err, marker="o", capsize=3, lw=1.8, color=color, label=label)
    ax.set_xticks(x, [DISPLAY_GENE[gene] for gene in TARGET_GENES])


def plot_main_figure(
    score_summary: pd.DataFrame,
    mutation_summary: pd.DataFrame,
    predicted_summary: pd.DataFrame,
    contrasts: pd.DataFrame,
    associations: pd.DataFrame,
    output: str | Path,
) -> None:
    fig, axes = plt.subplots(1, 5, figsize=(19, 4.4))
    ax = axes[0]
    ax.axis("off")
    ax.text(0.5, 0.88, "A. Mechanism", ha="center", fontweight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.62, "A$_H$>G$_H$ / C$_H$>T$_H$", ha="center", transform=ax.transAxes)
    ax.annotate("", xy=(0.5, 0.48), xytext=(0.5, 0.57), arrowprops={"arrowstyle": "->"}, xycoords="axes fraction")
    ax.text(0.5, 0.40, "real CDS codon", ha="center", transform=ax.transAxes)
    ax.annotate("", xy=(0.5, 0.26), xytext=(0.5, 0.35), arrowprops={"arrowstyle": "->"}, xycoords="axes fraction")
    ax.text(0.5, 0.18, "amino-acid change", ha="center", transform=ax.transAxes)

    ax = axes[1]
    colors = {"observed_A_rich_score": "#E45756", "observed_C_rich_score": "#4C78A8"}
    labels = {"observed_A_rich_score": "A-rich score", "observed_C_rich_score": "C-rich score"}
    for metric in colors:
        _errorbar(ax, score_summary.loc[score_summary["metric"].eq(metric)], "mean", labels[metric], colors[metric])
    ax.set_title("B. Observed composition")
    ax.set_ylabel("Mean codon-base richness score")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2]
    for mutation, color in zip(MUTATION_TYPES, ("#F58518", "#54A24B")):
        _errorbar(ax, mutation_summary.loc[mutation_summary["Mut"].eq(mutation)], "mean", MUTATION_DISPLAY[mutation], color)
    ax.set_title("C. Mutational pressure")
    ax.set_ylabel("Mean MutSpec")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[3]
    top_amino_acids = (
        predicted_summary.groupby("amino_acid")["mean"].agg(lambda values: np.nanmax(np.abs(values))).nlargest(5).index
    )
    palette = plt.get_cmap("tab10")
    for index, aa in enumerate(top_amino_acids):
        _errorbar(ax, predicted_summary.loc[predicted_summary["amino_acid"].eq(aa)], "mean", aa, palette(index))
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_title("D. Predicted AA pressure")
    ax.set_ylabel("Weighted normalized net opportunity")
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[4]
    label = "CytB - COX1"
    table = contrasts.loc[contrasts["contrast"].eq(label)]
    ax.axhline(0, color="0.8", lw=0.8)
    ax.axvline(0, color="0.8", lw=0.8)
    ax.scatter(table["predicted"], table["observed"], s=28, color="#7A5195")
    for row in table.itertuples(index=False):
        ax.annotate(row.amino_acid, (row.predicted, row.observed), xytext=(3, 2), textcoords="offset points", fontsize=7)
    association = associations.loc[associations["contrast"].eq(label)].iloc[0]
    ax.set_title("E. Predicted vs observed")
    ax.set_xlabel("Predicted change")
    ax.set_ylabel("Observed fraction change")
    ax.text(0.04, 0.96, f"Spearman ρ = {association.spearman_rho:.2f}\n95% CI [{association.ci_low:.2f}, {association.ci_high:.2f}]", transform=ax.transAxes, va="top", fontsize=8)

    fig.suptitle("COX1 → COX3 → CytB: observed composition and mutation-implied pressure", y=1.03, fontsize=13)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _write_csv(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, index=False)


__all__ = [
    "AA_ORDER", "BIOPYTHON_AVAILABLE", "DISPLAY_GENE", "MUTATION_DISPLAY",
    "TARGET_GENES", "EntrezSettings", "NCBITaxonomyClient",
    "amino_acid_codon_richness", "amino_acid_net_opportunities",
    "build_availability_and_attrition", "cluster_bootstrap_mean",
    "codon_mutational_opportunities", "load_mutation_cohort",
    "load_source_species_taxids", "mutation_weighted_shifts",
    "observed_amino_acid_composition", "observed_predicted_contrasts",
    "plot_main_figure", "safe_species_query",
    "sequence_qc", "translate_table2",
]
