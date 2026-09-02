# Signature Analysis: CO1, CO3, Cytb

Analysis of COSMIC v3.3 signature assignments in 52 matched mammalian species across three mitochondrial genes using SBS96 spectra.

## Workflow

Run both steps at once with `python run_all.py`, or one at a time as below.
Everything runs with the interpreter that starts the script, so activating this
folder's environment first is enough; `Rscript` must be on PATH.

### 1. Run SigProfilerAssignment
```bash
python 1_run_sigprofiler.py
```

This script:
- Prepares the SBS96 spectra from stage-5 matched data
- Runs COSMIC v3.3 assignment for low, high, and high-minus-low transition branches
- Aggregates results and extracts signature priors

**Output:** `results/sigprofiler_*.csv`

### 2. Run mSigAct Analysis
```bash
python 2_run_msigact.py
```

This script:
- Combines prepared branch matrices into a single mSigAct input
- Runs the R-based mSigAct analysis (if R is installed)
- Verifies results were created

**Output:** `results/msigact_*.csv`

### 3. Generate Report and Figures
Open and run the Jupyter notebook:
```bash
jupyter notebook SignatureAnalysis.ipynb
```

This notebook:
- Loads results from both SigProfilerAssignment and mSigAct
- Normalizes activities and calculates cross-tool concordance
- Generates publication-style figures comparing three genes
- Produces summary statistics

**Outputs:**
- `results/tool_activity_fractions.csv` — normalized activities
- `results/tool_concordance.csv` — cross-tool agreement metrics
- `figures/signature_activities_by_gene.png` — main comparison figure
- `figures/signature_activities_by_gene.pdf` — publication-ready PDF

## Directory Structure

```
signatures/
├── run_all.py                    # Cross-platform runner for steps 1 and 2
├── 1_run_sigprofiler.py          # SigProfilerAssignment runner
├── 2_run_msigact.py              # mSigAct runner
├── run_msigact.R                 # R script for mSigAct (called by 2_run_msigact.py)
├── signature_analysis.py         # Shared analysis utilities
├── SignatureAnalysis.ipynb       # Main analysis notebook
│
├── reference/                    # Human triplet counts used to renormalize
│   └── triplet_counts_GRCh37_upper.json
│
├── input/                        # Prepared SBS96 matrices (rebuilt by step 1)
│   ├── low_Ts_samples.txt
│   ├── high_Ts_samples.txt
│   ├── high_minus_low_Ts_samples.txt
│   ├── samples_mSigAct.txt
│   └── input_manifest.csv
│
├── results/                      # Final analysis results
│   ├── sigprofiler_activities.csv
│   ├── sigprofiler_priors.csv
│   ├── sigprofiler_quality.csv
│   ├── msigact_activities.csv
│   ├── msigact_quality.csv
│   ├── tool_activity_fractions.csv
│   ├── tool_concordance.csv
│   └── *.json / *.csv (metadata)
│
├── figures/                      # Generated figures
│   ├── signature_activities_by_gene.png
│   └── signature_activities_by_gene.pdf
│
└── work/                         # Intermediate files (not tracked)
    ├── sigprofiler/
    └── msigact/
```

## Requirements

### Python

This substage pins its own environment because SigProfilerAssignment constrains
NumPy and pandas differently from the rest of the repository:

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv/Scripts/activate
python -m pip install -r requirements.txt
```

### R (optional, for mSigAct)

mSigAct is an R package. If you want to run mSigAct analysis, install R and run:
```bash
Rscript install_r_dependencies.R
```

Or manually in R:
```r
install.packages(c("ICAMS", "cosmicsig", "mSigAct"))
```

## Configuration

Both assignment tools analyze three genes with two spectrum types:

- **Genes:** CO1, CO3, Cytb
- **Branches:** low_Ts, high_Ts, high_minus_low_Ts
- **Spectrum types:** Ts_only (transitions only), Ts_and_Tv (transitions + transversions)
- **Total samples:** 3 × 3 × 2 = 18 per tool

## Notes

- Input data is derived from stage-5 matched 192-component spectra
- COSMIC v3.3 (GRCh37 context renormalization)
- Cross-tool concordance measured via cosine similarity and total variation distance
- `input/` and `work/` are rebuilt by the scripts and are not version-controlled;
  `reference/`, `results/`, and `figures/` are tracked so a fresh clone can
  check the reported numbers without re-running the tools
