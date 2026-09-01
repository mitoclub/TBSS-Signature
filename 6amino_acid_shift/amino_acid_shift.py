"""Minimal sequence-QC helpers for the MIDORI2 amino-acid ratio analysis.

All statistical analysis intentionally lives in ``AminoAcidShift.ipynb``.
This module contains only constants, vertebrate mitochondrial translation,
and validation of the paired nucleotide/protein records supplied by MIDORI2.
"""

from __future__ import annotations

import pandas as pd

from .utils import TARGET_GENES

# Broad mammalian protein-length screens. MIDORI calls ATP6/ATP8 A6/A8.
PROTEIN_LENGTH_BOUNDS = {
    "CO1": (450, 550),
    "CO2": (200, 300),
    "A8": (35, 100),
    "A6": (180, 280),
    "CO3": (220, 300),
    "ND3": (90, 140),
    "ND4L": (70, 120),
    "ND4": (400, 500),
    "ND5": (550, 650),
    "Cytb": (330, 420),
}

# Vertebrate mitochondrial genetic code (NCBI translation table 2).
TABLE2 = {
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
TABLE2_STARTS = frozenset(("ATT", "ATC", "ATA", "ATG", "GTG"))


def translate_table2(coding_sequence: str, *, apply_start_rule: bool = True) -> str:
    """Translate complete codons with vertebrate mitochondrial table 2."""

    sequence = str(coding_sequence).upper()
    codons = [sequence[index:index + 3] for index in range(0, len(sequence) - 2, 3)]
    protein = [TABLE2.get(codon, "X") for codon in codons]
    if apply_start_rule and codons and codons[0] in TABLE2_STARTS:
        protein[0] = "M"
    return "".join(protein)


def sequence_qc(
    sequence_data: pd.DataFrame,
    *,
    maximum_ambiguous_fraction: float = 0.01,
) -> pd.DataFrame:
    """Validate selected MIDORI2 CDS/protein pairs under translation table 2."""

    if sequence_data.empty:
        columns = list(sequence_data.columns) + [
            "coding_cds_sequence",
            "computed_translation",
            "translation_match",
            "terminal_incomplete_stop",
            "sequence_qc_status",
            "sequence_qc_reason",
        ]
        return pd.DataFrame(columns=list(dict.fromkeys(columns)))

    rows: list[tuple[dict[str, object], list[str]]] = []
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
        if bool(row.get("is_partial")):
            reasons.append("partial_CDS")

        ambiguous = sum(base not in "ACGT" for base in coding)
        if coding and ambiguous / len(coding) > maximum_ambiguous_fraction:
            reasons.append("excess_ambiguous_nucleotides")
        remainder = len(coding) % 3
        terminal_fragment = coding[-remainder:] if remainder else ""
        terminal_incomplete = terminal_fragment in {"T", "TA"}
        if remainder and not terminal_incomplete:
            reasons.append("frameshift_or_unrecognized_terminal_fragment")

        calculated = translate_table2(
            coding,
            apply_start_rule=bool(row.get("apply_start_rule", True)),
        )
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
    for index, (row, reasons) in enumerate(rows):
        low, high = PROTEIN_LENGTH_BOUNDS[str(row["Gene"])]
        protein_length = float(row.get("protein_length", 0))
        if not (low <= protein_length <= high):
            reasons.append(f"protein_length_outlier[{low},{high}]")
        result.loc[index, "sequence_qc_status"] = "pass" if not reasons else "fail"
        result.loc[index, "sequence_qc_reason"] = "pass" if not reasons else ";".join(reasons)
    return result


__all__ = [
    "TABLE2",
    "PROTEIN_LENGTH_BOUNDS",
    "TARGET_GENES",
    "sequence_qc",
    "translate_table2",
]
