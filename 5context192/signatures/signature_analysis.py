"""Prepare and aggregate COSMIC SBS96 signature assignments by gene.

The stage-5 source contains strand-specific SBS192 spectra.  COSMIC v3.3
signatures used by SigProfilerAssignment and mSigAct are SBS96 profiles, so the
conversion performed here is deliberate and strand aware.  It follows the
``4signatures`` workflow in the reference repository:

* average complete profiles over the fixed matched-species cohort;
* separate the high (C>T, A>G) and low (G>A, T>C) heavy-strand transitions;
* reverse-complement non-C/T-centered components while collapsing to SBS96;
* average the two members of every transversion pair;
* create high, low, and positive high-minus-low spectra, with and without
  transversions; and
* rescale the profiles with the same uppercase GRCh37 trinucleotide counts and
  coefficient used by the reference analysis.

No signature tool is imported by this module, which keeps preparation and QC
testable without either heavy optional dependency.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


GENES = ("CO1", "CO3", "Cytb")
GENE_LABELS = {"CO1": "COX1", "CO3": "COX3", "Cytb": "CytB"}
BRANCHES = ("low_Ts", "high_Ts", "high_minus_low_Ts")
BRANCH_LABELS = {
    "low_Ts": "Low transitions",
    "high_Ts": "High transitions",
    "high_minus_low_Ts": "High - low transitions",
}
SPECTRUM_TYPES = ("Ts_only", "Ts_and_Tv")

BASES = ("A", "C", "G", "T")
COSMIC_SBS6 = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
SBS96_ORDER = tuple(
    f"{left}[{substitution}]{right}"
    for substitution in COSMIC_SBS6
    for left in BASES
    for right in BASES
)
HIGH_TRANSITIONS = frozenset(("C>T", "A>G"))
LOW_TRANSITIONS = frozenset(("G>A", "T>C"))
TRANSVERSIONS = frozenset(
    ("A>C", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>G")
)
COMPLEMENT = str.maketrans("ACGT", "TGCA")
HUMAN_SCALE_COEFFICIENT = 6.6e-5


@dataclass(frozen=True)
class ParsedContext:
    left: str
    reference: str
    alternate: str
    right: str

    @property
    def substitution(self) -> str:
        return f"{self.reference}>{self.alternate}"

    @property
    def trinucleotide(self) -> str:
        return f"{self.left}{self.reference}{self.right}"

    def __str__(self) -> str:
        return f"{self.left}[{self.reference}>{self.alternate}]{self.right}"


def parse_context(label: str) -> ParsedContext:
    """Parse one exact ``L[REF>ALT]R`` label without a regex dependency."""

    text = str(label)
    if len(text) != 7 or text[1] != "[" or text[3] != ">" or text[5] != "]":
        raise ValueError(f"Invalid context label: {label!r}")
    parsed = ParsedContext(text[0], text[2], text[4], text[6])
    if any(base not in BASES for base in (
        parsed.left, parsed.reference, parsed.alternate, parsed.right
    )):
        raise ValueError(f"Invalid context label: {label!r}")
    if parsed.reference == parsed.alternate:
        raise ValueError(f"Reference and alternate must differ: {label!r}")
    return parsed


def reverse_complement_context(label: str) -> str:
    parsed = parse_context(label)
    return str(
        ParsedContext(
            parsed.right.translate(COMPLEMENT),
            parsed.reference.translate(COMPLEMENT),
            parsed.alternate.translate(COMPLEMENT),
            parsed.left.translate(COMPLEMENT),
        )
    )


def canonical_sbs96_context(label: str) -> str:
    """Return the COSMIC C/T-centred representative of an SBS192 context."""

    parsed = parse_context(label)
    canonical = label if parsed.reference in "CT" else reverse_complement_context(label)
    if canonical not in SBS96_ORDER:
        raise ValueError(f"Could not map {label!r} to a canonical SBS96 context")
    return canonical


def load_matched_heavy_spectra(
    path: str | Path,
    *,
    genes: Sequence[str] = GENES,
    expected_class: str = "Mammalia",
) -> pd.DataFrame:
    """Load and validate the complete, already oriented stage-5 cohort."""

    data = pd.read_csv(path)
    required = {"Gene", "Class", "Species", "Mut", "MutSpec", "ContextOrientation"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Matched SBS192 table is missing columns: {sorted(missing)}")
    data = data.loc[data["Gene"].isin(genes)].copy()
    if data.empty:
        raise ValueError("Matched SBS192 table contains no requested genes")
    if set(data["Gene"]) != set(genes):
        absent = sorted(set(genes).difference(data["Gene"]))
        raise ValueError(f"Matched SBS192 table is missing genes: {absent}")
    if set(data["Class"]) != {expected_class}:
        raise ValueError(
            f"Expected only {expected_class}; observed {sorted(set(data['Class']))}"
        )
    if set(data["ContextOrientation"].astype(str).str.lower()) != {"heavy"}:
        raise ValueError("Signature input must already be heavy-strand oriented")

    values = pd.to_numeric(data["MutSpec"], errors="coerce")
    if values.isna().any() or not np.isfinite(values).all() or values.lt(0).any():
        raise ValueError("MutSpec must be finite and non-negative")
    data["MutSpec"] = values
    duplicated = data.duplicated(["Species", "Gene", "Mut"])
    if duplicated.any():
        raise ValueError("Duplicate Species/Gene/Mut rows in matched SBS192 table")

    profile_sizes = data.groupby(["Species", "Gene"], sort=False).size()
    if not profile_sizes.eq(192).all():
        raise ValueError("Every matched Species/Gene profile must contain 192 rows")
    profile_sums = data.groupby(["Species", "Gene"], sort=False)["MutSpec"].sum()
    if not np.allclose(profile_sums.to_numpy(), 1.0, atol=1e-8, rtol=0):
        raise ValueError("Every matched Species/Gene MutSpec profile must sum to one")
    per_species = data.groupby("Species", sort=False)["Gene"].nunique()
    if not per_species.eq(len(genes)).all():
        raise ValueError("The signature cohort must contain every requested gene per species")
    return data


def mean_gene_sbs192(
    data: pd.DataFrame,
    *,
    genes: Sequence[str] = GENES,
) -> pd.DataFrame:
    """Calculate equally species-weighted mean SBS192 profiles by gene."""

    mean_long = (
        data.groupby(["Gene", "Mut"], sort=False, observed=True)["MutSpec"]
        .mean()
        .rename("Weight")
        .reset_index()
    )
    sizes = mean_long.groupby("Gene", sort=False).size().reindex(genes)
    if sizes.isna().any() or not sizes.eq(192).all():
        raise ValueError("Each mean gene spectrum must contain 192 contexts")
    sums = mean_long.groupby("Gene", sort=False)["Weight"].sum().reindex(genes)
    if not np.allclose(sums.to_numpy(), 1.0, atol=1e-8, rtol=0):
        raise ValueError("Mean gene SBS192 profiles must sum to one")
    return mean_long


def collapse_selected_components(
    mean_long: pd.DataFrame,
    substitutions: frozenset[str],
    *,
    genes: Sequence[str] = GENES,
) -> pd.DataFrame:
    """Select SBS192 substitutions and collapse complementary contexts."""

    selected_rows: list[dict[str, object]] = []
    for row in mean_long.itertuples(index=False):
        parsed = parse_context(row.Mut)
        if parsed.substitution not in substitutions:
            continue
        selected_rows.append(
            {
                "Gene": row.Gene,
                "MutationType": canonical_sbs96_context(row.Mut),
                "Weight": float(row.Weight),
            }
        )
    if not selected_rows:
        raise ValueError(f"No components found for substitutions={sorted(substitutions)}")
    collapsed = (
        pd.DataFrame(selected_rows)
        .groupby(["Gene", "MutationType"], observed=True)["Weight"]
        .sum()
        .unstack(fill_value=0.0)
        .reindex(index=list(genes), columns=list(SBS96_ORDER), fill_value=0.0)
    )
    collapsed.index.name = "Gene"
    collapsed.columns.name = "MutationType"
    return collapsed.astype(float)


def build_spectrum_variants(
    mean_long: pd.DataFrame,
    *,
    genes: Sequence[str] = GENES,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    """Build the low, high, and positive high-minus-low SBS96 variants."""

    low = collapse_selected_components(mean_long, LOW_TRANSITIONS, genes=genes)
    high = collapse_selected_components(mean_long, HIGH_TRANSITIONS, genes=genes)
    transversions = (
        collapse_selected_components(mean_long, TRANSVERSIONS, genes=genes) / 2.0
    )
    difference = (high - low).clip(lower=0.0)
    variants = {
        "low_Ts": (low, low + transversions),
        "high_Ts": (high, high + transversions),
        "high_minus_low_Ts": (difference, difference + transversions),
    }
    for branch, pair in variants.items():
        for spectrum in pair:
            if spectrum.shape != (len(genes), 96):
                raise AssertionError(f"Unexpected {branch} matrix shape: {spectrum.shape}")
            if not np.isfinite(spectrum.to_numpy()).all() or spectrum.lt(0).any().any():
                raise ValueError(f"Invalid values in {branch}")
    return variants


def load_human_triplet_counts(path: str | Path) -> dict[str, int]:
    counts = json.loads(Path(path).read_text(encoding="utf-8"))
    expected = {f"{a}{b}{c}" for a in BASES for b in BASES for c in BASES}
    if set(counts) != expected:
        missing = sorted(expected.difference(counts))
        extra = sorted(set(counts).difference(expected))
        raise ValueError(f"GRCh37 triplet counts mismatch; missing={missing}, extra={extra}")
    parsed = {key: int(value) for key, value in counts.items()}
    if any(value <= 0 for value in parsed.values()):
        raise ValueError("All GRCh37 triplet counts must be positive integers")
    return parsed


def rescale_to_human_counts(
    spectra: pd.DataFrame,
    triplet_counts: Mapping[str, int],
    *,
    scale_coefficient: float = HUMAN_SCALE_COEFFICIENT,
) -> pd.DataFrame:
    """Turn profile weights into integer human-context-emulated counts."""

    if list(spectra.columns) != list(SBS96_ORDER):
        raise ValueError("SBS96 columns are missing or not in canonical order")
    multipliers = np.array(
        [triplet_counts[parse_context(label).trinucleotide] for label in SBS96_ORDER],
        dtype=float,
    )
    scaled = np.rint(
        spectra.to_numpy(dtype=float) * multipliers[None, :] * scale_coefficient
    ).astype(np.int64)
    if (scaled < 0).any() or scaled.sum(axis=1).min() <= 0:
        raise ValueError("Rescaled spectra must have positive non-negative integer counts")
    return pd.DataFrame(
        scaled,
        index=spectra.index.copy(),
        columns=spectra.columns.copy(),
    )


def prepare_signature_inputs(
    matched_path: str | Path,
    triplet_counts_path: str | Path,
    output_dir: str | Path,
    *,
    genes: Sequence[str] = GENES,
) -> pd.DataFrame:
    """Create three tab-delimited SigProfiler/mSigAct input matrices."""

    genes = tuple(genes)
    matched = load_matched_heavy_spectra(matched_path, genes=genes)
    mean_long = mean_gene_sbs192(matched, genes=genes)
    variants = build_spectrum_variants(mean_long, genes=genes)
    triplet_counts = load_human_triplet_counts(triplet_counts_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    for branch in BRANCHES:
        ts_only, ts_and_tv = variants[branch]
        sample_rows = []
        sample_names = []
        for spectrum_type, spectrum in zip(SPECTRUM_TYPES, (ts_only, ts_and_tv)):
            scaled = rescale_to_human_counts(spectrum, triplet_counts)
            for gene in genes:
                sample_name = f"{gene}__{spectrum_type}"
                sample_names.append(sample_name)
                sample_rows.append(scaled.loc[gene].to_numpy())
                manifest_rows.append(
                    {
                        "Branch": branch,
                        "Sample": sample_name,
                        "Gene": gene,
                        "GeneLabel": GENE_LABELS.get(gene, gene),
                        "SpectrumType": spectrum_type,
                        "IncludesTransversions": spectrum_type == "Ts_and_Tv",
                        "NMatchedSpecies": int(matched["Species"].nunique()),
                        "TotalIntegerCounts": int(scaled.loc[gene].sum()),
                    }
                )
        matrix = pd.DataFrame(sample_rows, index=sample_names, columns=SBS96_ORDER).T
        matrix.index.name = "MutationType"
        matrix.to_csv(output / f"{branch}_samples.txt", sep="\t")

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(output / "input_manifest.csv", index=False)
    return manifest


def parse_sample_name(sample: str) -> tuple[str, str]:
    parts = str(sample).split("__")
    if len(parts) != 2 or parts[1] not in SPECTRUM_TYPES:
        raise ValueError(f"Unexpected sample name: {sample!r}")
    return parts[0], parts[1]


def parse_msigact_sample_name(sample: str) -> tuple[str, str, str]:
    parts = str(sample).split("__")
    if len(parts) != 3 or parts[0] not in BRANCHES or parts[2] not in SPECTRUM_TYPES:
        raise ValueError(f"Unexpected mSigAct sample name: {sample!r}")
    return parts[0], parts[1], parts[2]


__all__ = [
    "BRANCHES",
    "BRANCH_LABELS",
    "COSMIC_SBS6",
    "GENES",
    "GENE_LABELS",
    "HUMAN_SCALE_COEFFICIENT",
    "SBS96_ORDER",
    "SPECTRUM_TYPES",
    "build_spectrum_variants",
    "canonical_sbs96_context",
    "collapse_selected_components",
    "load_human_triplet_counts",
    "load_matched_heavy_spectra",
    "mean_gene_sbs192",
    "parse_context",
    "parse_msigact_sample_name",
    "parse_sample_name",
    "prepare_signature_inputs",
    "rescale_to_human_counts",
    "reverse_complement_context",
]
