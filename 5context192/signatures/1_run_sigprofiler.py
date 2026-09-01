"""Run COSMIC v3.3 assignment with SigProfilerAssignment for all gene spectra.

Usage:
  python 1_run_sigprofiler.py
"""

from pathlib import Path
import json
import pandas as pd

# Import from signature_analysis module
from signature_analysis import BRANCHES, GENES, GENE_LABELS, SPECTRUM_TYPES, prepare_signature_inputs, parse_sample_name

HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "input"
WORK_DIR = HERE / "work" / "sigprofiler"
RESULTS_DIR = HERE / "results"

EXCLUDED_SIGNATURE_SUBGROUPS = [
    "Artifact_signatures",
    "Immunosuppressants_signatures",
    "Treatment_signatures",
    "Lymphoid_signatures",
    "Colibactin_signatures",
    "AA_signatures",
]


def prepare_inputs():
    """Prepare SBS96 spectra from stage-5 matched data."""
    stage = HERE.parent
    manifest = prepare_signature_inputs(
        stage / "data" / "matched_192_spectra_CO1_CO3_Cytb.csv.gz",
        HERE / "reference" / "triplet_counts_GRCh37_upper.json",
        INPUT_DIR,
        genes=GENES,
    )
    print(
        f"Prepared {len(manifest)} gene/branch/spectrum inputs from "
        f"{manifest['NMatchedSpecies'].iloc[0]} matched mammalian species."
    )
    return manifest


def run_assignments(*, make_plots: bool = True, cpu: int = 1) -> str:
    """Run SigProfilerAssignment for all prepared gene spectra."""
    try:
        import SigProfilerAssignment
        from SigProfilerAssignment import Analyzer as Analyze
    except ImportError as error:
        raise RuntimeError(
            "SigProfilerAssignment is not installed. "
            "Install it with: pip install SigProfilerAssignment"
        ) from error

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    for branch in BRANCHES:
        samples = INPUT_DIR / f"{branch}_samples.txt"
        if not samples.is_file():
            raise FileNotFoundError(
                f"Missing {samples}; run prepare_inputs first"
            )
        output = WORK_DIR / branch
        output.mkdir(parents=True, exist_ok=True)
        print(f"SigProfilerAssignment: {branch}")
        Analyze.cosmic_fit(
            samples=str(samples),
            output=str(output),
            input_type="matrix",
            context_type="96",
            genome_build="GRCh37",
            cosmic_version=3.3,
            exome=False,
            collapse_to_SBS96=True,
            nnls_remove_penalty=0.01,
            nnls_add_penalty=0.02,
            exclude_signature_subgroups=EXCLUDED_SIGNATURE_SUBGROUPS,
            export_probabilities=False,
            make_plots=make_plots,
            sample_reconstruction_plots=False,
            verbose=False,
            cpu=cpu,
            add_background_signatures=True,
        )
    return str(SigProfilerAssignment.__version__)


def aggregate_assignments() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Aggregate SigProfilerAssignment results from all branches."""
    activity_frames = []
    quality_frames = []

    for branch in BRANCHES:
        solution = WORK_DIR / branch / "Assignment_Solution"
        activities_path = solution / "Activities" / "Assignment_Solution_Activities.txt"
        stats_path = solution / "Solution_Stats" / "Assignment_Solution_Samples_Stats.txt"

        if not activities_path.is_file() or not stats_path.is_file():
            raise FileNotFoundError(f"Incomplete results for {branch}: {solution}")

        # Parse activities
        wide = pd.read_csv(activities_path, sep="\t", index_col=0)
        wide.index = wide.index.astype(str)
        wide.index.name = "Sample"
        long = wide.reset_index().melt(
            id_vars="Sample", var_name="Signature", value_name="Activity"
        )
        long["Activity"] = pd.to_numeric(long["Activity"], errors="raise")
        long["Branch"] = branch
        parsed = long["Sample"].map(parse_sample_name)
        long["Gene"] = parsed.str[0]
        long["SpectrumType"] = parsed.str[1]
        denominator = long.groupby(["Branch", "Sample"])["Activity"].transform("sum")
        long["ActivityFraction"] = long["Activity"].div(
            denominator.where(denominator.gt(0))
        ).fillna(0.0)
        activity_frames.append(long)

        # Parse quality stats
        stats = pd.read_csv(stats_path, sep="\t")
        sample_column = "Sample Names" if "Sample Names" in stats.columns else stats.columns[0]
        stats = stats.rename(columns={sample_column: "Sample"})
        stats["Sample"] = stats["Sample"].astype(str)
        stats["Branch"] = branch
        parsed_stats = stats["Sample"].map(parse_sample_name)
        stats["Gene"] = parsed_stats.str[0]
        stats["SpectrumType"] = parsed_stats.str[1]
        quality_frames.append(stats)

    activities = pd.concat(activity_frames, ignore_index=True)
    activities = activities.loc[activities["Activity"].gt(0)].reset_index(drop=True)
    quality = pd.concat(quality_frames, ignore_index=True)

    # Calculate priors
    signature_totals = activities.groupby("Signature", as_index=False)["Activity"].sum()
    total_activity = signature_totals["Activity"].sum()
    if total_activity <= 0:
        raise ValueError("SigProfilerAssignment returned no assigned activity")

    signature_totals["PriorUnrounded"] = signature_totals["Activity"] / total_activity
    signature_totals["Prior"] = signature_totals["PriorUnrounded"].round(2)
    priors = (
        signature_totals.loc[
            signature_totals["Activity"].gt(0) & signature_totals["Prior"].ge(0.01)
        ]
        .rename(columns={"Activity": "TotalActivity"})
        .sort_values(["Prior", "Signature"], ascending=[False, True])
        .reset_index(drop=True)
    )

    if priors.empty:
        raise ValueError("No SigProfiler signature reached 1% prior threshold")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    activities.to_csv(RESULTS_DIR / "sigprofiler_activities.csv", index=False)
    quality.to_csv(RESULTS_DIR / "sigprofiler_quality.csv", index=False)
    priors.to_csv(RESULTS_DIR / "sigprofiler_priors.csv", index=False)

    return activities, quality, priors


def write_metadata(version: str, *, make_plots: bool, cpu: int) -> None:
    """Write metadata about the SigProfilerAssignment run."""
    metadata = {
        "tool": "SigProfilerAssignment",
        "version": version,
        "cosmic_version": "3.3",
        "genome_build": "GRCh37",
        "input_type": "matrix",
        "context_type": "96",
        "source_context_type": "192",
        "nnls_remove_penalty": 0.01,
        "nnls_add_penalty": 0.02,
        "add_background_signatures": True,
        "excluded_signature_subgroups": EXCLUDED_SIGNATURE_SUBGROUPS,
        "make_plots": make_plots,
        "cpu": cpu,
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "sigprofiler_run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    print("Step 1: Preparing inputs...")
    manifest = prepare_inputs()
    print(f"  ✓ {len(manifest)} spectra prepared")

    print("\nStep 2: Running SigProfilerAssignment...")
    version = run_assignments(make_plots=True, cpu=-1)
    print(f"  ✓ SigProfilerAssignment {version}")

    print("\nStep 3: Aggregating results...")
    activities, quality, priors = aggregate_assignments()
    write_metadata(version, make_plots=True, cpu=-1)
    print(f"  ✓ Aggregated {activities[['Branch', 'Sample']].drop_duplicates().shape[0]} spectra")
    print(f"  ✓ Found {len(priors)} signatures (≥1% prior)")
    print("\nTop signatures:")
    print(priors[["Signature", "Prior", "TotalActivity"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
