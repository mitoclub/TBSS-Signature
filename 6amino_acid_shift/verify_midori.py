"""Verify the minimal MIDORI2 data contract and notebook figures."""

from __future__ import annotations

from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from .utils import TARGET_GENES, sha256


PSEUDOCOUNT = 0.5
OBSOLETE_OUTPUTS = {
    "amino_acid_counts_by_species_gene.csv",
    "amino_acid_ratio_monotonic_trend.csv",
    "amino_acid_ratio_paired_tests.csv",
    "amino_acid_ratio_summary_by_gene.csv",
    "amino_acid_ratios_by_species_gene.csv",
    "analysis_species.csv",
    "midori_record_consistency.csv",
    "midori_sequence_manifest.csv.gz",
    "midori_species_catalog.csv",
    "sequence_qc.csv.gz",
    "species_gene_availability.csv",
    "species_inclusion_status.csv",
    "species_intersection_summary.csv",
}


def _ratios(row: object) -> tuple[float, float]:
    amino_acids = Counter(str(row.computed_translation).rstrip("*"))
    coding = str(row.coding_cds_sequence)
    codons = [coding[index:index + 3] for index in range(0, len(coding) - 2, 3)]
    leucine_ttr = sum(codon in {"TTA", "TTG"} for codon in codons)
    return (
        (amino_acids["N"] + amino_acids["K"] + PSEUDOCOUNT)
        / (amino_acids["G"] + PSEUDOCOUNT),
        (amino_acids["P"] + PSEUDOCOUNT)
        / (amino_acids["F"] + leucine_ttr + PSEUDOCOUNT),
    )


def main() -> None:
    stage = Path(__file__).resolve().parent
    derived = stage / "data" / "derived_midori"
    cache = stage / "data" / "midori"

    files = pd.read_csv(derived / "midori_file_manifest.csv", dtype={"sha256": str})
    assert len(files) == 2 * len(TARGET_GENES)
    assert set(files["Gene"]) == TARGET_GENES
    assert set(files["molecule"]) == {"NUC", "AA"}
    assert files[["Gene", "molecule"]].duplicated().sum() == 0
    for row in files.itertuples(index=False):
        path = cache / row.filename
        assert path.exists() and path.stat().st_size == int(row.compressed_bytes)
        assert sha256(path) == row.sha256
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            assert handle.readline().startswith(">")

    sequences = pd.read_csv(
        derived / "midori_gene_sequences.csv.gz",
        dtype={"TaxID": str},
    )
    required = {
        "Species", "TaxID", "Gene", "coding_cds_sequence",
        "computed_translation", "sequence_qc_status", "translation_match",
        "same_feature_pair", "transl_table",
    }
    assert required.issubset(sequences.columns)
    assert not sequences.duplicated(["Species", "Gene"]).any()
    assert set(sequences["Gene"]) == TARGET_GENES
    assert sequences.groupby("Species")["Gene"].nunique().eq(len(TARGET_GENES)).all()
    assert sequences.groupby("Species")["TaxID"].nunique().eq(1).all()
    assert sequences["sequence_qc_status"].eq("pass").all()
    assert sequences["translation_match"].eq(True).all()
    assert sequences["same_feature_pair"].eq(True).all()
    assert sequences["transl_table"].eq(2).all()

    ratio_rows = []
    for row in sequences.itertuples(index=False):
        ratio_rows.append((row.Species, row.Gene, *_ratios(row)))
    ratios = pd.DataFrame(
        ratio_rows,
        columns=["Species", "Gene", "(Asn+Lys)/Gly", "Pro/(Phe+LeuTTR)"],
    )
    assert np.isfinite(ratios.iloc[:, 2:].to_numpy()).all()
    for ratio_name in ratios.columns[2:]:
        wide = ratios.pivot(index="Species", columns="Gene", values=ratio_name)
        assert np.median(wide["Cytb"] / wide["CO1"]) > 1

    notebook_path = stage / "AminoAcidShift.ipynb"
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    assert notebook["nbformat"] == 4
    code = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            compile("".join(cell["source"]), f"notebook-cell-{index}", "exec")
    assert ".to_csv(" not in code
    assert "friedmanchisquare" in code
    assert "page_trend_test" in code
    assert "Kendall's W" in code
    assert "holm_adjust" in code
    assert "p-value Holm (45 pairs)" in code
    assert "wilcoxon" in code.casefold()
    assert "alternative='greater'" in code
    assert "alternative='two-sided'" not in code
    assert "combinations(range(len(GENE_ORDER)), 2)" in code
    assert "counts['Pro/Phe']" not in code
    assert "PSEUDOCOUNT = 0.5" in code
    assert "paired wilcoxon, one-sided" in code.casefold()
    assert "+ 'p-value ' + format_p_value(endpoint_p)" in code
    assert "$H_1$" not in code
    assert "bootstrap" not in code.casefold()

    for stem in ("amino_acid_ratios_by_gene", "paired_ratio_fold_changes"):
        png = stage / "figures" / f"{stem}.png"
        pdf = stage / "figures" / f"{stem}.pdf"
        assert png.exists() and png.stat().st_size > 50_000
        assert pdf.exists() and pdf.stat().st_size > 10_000
        with Image.open(png) as image:
            assert image.width >= 2_000 and image.height >= 1_000
        assert pdf.read_bytes()[:4] == b"%PDF"

    assert not any((derived / name).exists() for name in OBSOLETE_OUTPUTS)
    assert not (stage / "figures" / "midori" / "main_amino_acid_shift.png").exists()

    print(
        "MIDORI2 verification passed: "
        f"{sequences['Species'].nunique():,} paired species, "
        f"{len(TARGET_GENES)} Major Arc genes, two minimal data outputs, "
        "and two publication-style figures."
    )


if __name__ == "__main__":
    main()
