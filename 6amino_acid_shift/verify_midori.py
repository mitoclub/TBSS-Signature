"""Verify persisted invariants of the stage-6 MIDORI2 analysis."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd


def repository_root() -> Path:
    here = Path.cwd().resolve()
    for candidate in (here, *here.parents):
        if (candidate / "6amino_acid_shift" / "midori_analysis.py").exists():
            return candidate
    raise FileNotFoundError("Run from the repository root or 6amino_acid_shift")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    root = repository_root()
    stage = root / "6amino_acid_shift"
    cache = stage / "data" / "midori"
    derived = stage / "data" / "derived_midori"

    files = pd.read_csv(derived / "midori_file_manifest.csv", dtype=str)
    assert len(files) == 6
    for row in files.itertuples(index=False):
        path = cache / row.filename
        assert path.exists(), path
        assert path.stat().st_size == int(row.compressed_bytes), path
        assert sha256(path) == row.sha256, path
        with gzip.open(path, "rb") as handle:
            assert handle.read(1) == b">", path

    availability = pd.read_csv(derived / "species_gene_availability.csv")
    qc = pd.read_csv(derived / "sequence_qc.csv")
    final = pd.read_csv(derived / "final_common_species.csv")
    matching = pd.read_csv(derived / "species_midori_matching.csv")
    consistency = pd.read_csv(derived / "midori_record_consistency.csv")
    manifest = pd.read_csv(derived / "midori_sequence_manifest.csv")
    canonical = pd.read_csv(derived / "midori_gene_sequences.csv.gz")

    assert len(availability) == 52
    assert availability[["has_seq_CO1", "has_seq_CO3", "has_seq_Cytb"]].all().all()
    assert int(availability["included_final"].sum()) == 51
    assert len(final) == 51 and final["Species"].is_unique
    assert set(final["Species"]) == set(
        availability.loc[availability["included_final"], "Species"]
    )
    assert len(qc) == len(canonical) == 52 * 3
    assert int(qc["sequence_qc_status"].eq("pass").sum()) == 51 * 3 + 2
    assert qc["same_feature_pair"].all()
    assert qc["translation_match"].all()
    final_qc = qc.loc[qc["Species"].isin(final["Species"])]
    assert len(final_qc) == 51 * 3
    assert final_qc["sequence_qc_status"].eq("pass").all()
    assert matching["match_status"].isin(["exact", "taxid_confirmed"]).all()
    assert int(consistency["same_accession_all_three"].sum()) == 50
    assert manifest.loc[manifest["selected"], ["Species", "Gene"]].duplicated().sum() == 0
    assert len(manifest.loc[manifest["selected"]]) == 52 * 3

    composition = pd.read_csv(derived / "observed_amino_acid_composition.csv")
    opportunities = pd.read_csv(derived / "codon_mutational_opportunities.csv.gz")
    predicted = pd.read_csv(derived / "mutation_weighted_amino_acid_shifts.csv")
    assert len(composition) == 51 * 3 * 20
    assert set(composition["Species"]) == set(final["Species"])
    assert set(opportunities["Species"]) == set(final["Species"])
    assert len(predicted) == 51 * 3 * 20
    assert predicted["combined_predicted_shift"].notna().all()

    notebook = json.loads((stage / "AminoAcidShift.ipynb").read_text(encoding="utf-8"))
    assert notebook["nbformat"] == 4
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            compile("".join(cell["source"]), f"notebook-cell-{index}", "exec")

    figure = stage / "figures" / "midori" / "main_amino_acid_shift.png"
    assert figure.exists() and figure.stat().st_size > 10_000, figure

    print("MIDORI2 verification passed: 52 species found, 51 final, 50 shared-accession.")


if __name__ == "__main__":
    main()
