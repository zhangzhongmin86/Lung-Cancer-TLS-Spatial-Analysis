# Lung-Cancer-TLS-Spatial-Analysis

Analysis code for the accompanying manuscript.

## Repository contents

- `R/main_figures/`: R scripts used for main figures 1 and 2.
- `R/supplementary_figures/`: R scripts used for supplementary figures.
- `notebooks/main_figures/`: Python notebooks used for main figures 4–6.
- `notebooks/supplementary_figures/`: Python notebooks used for supplementary analyses.
- `commune/`: scripts and notebooks implementing the COMMUNE workflow.

```text
.
├── R/
│   ├── main_figures/
│   └── supplementary_figures/
├── notebooks/
│   ├── main_figures/
│   └── supplementary_figures/
├── commune/
├── .gitignore
└── README.md
```

### Main analysis index

| Analysis | Code |
| --- | --- |
| Atlas construction and overview | `R/main_figures/figure1.R` |
| Cell-composition analysis | `R/main_figures/figure2.R` |
| TLS identification | `notebooks/main_figures/figure4-Identifying TLS.ipynb` |
| Tumor-region classification | `notebooks/main_figures/figure4-Regional grouping.ipynb` |
| TLS trajectory classification | `notebooks/main_figures/Fig5_TLS_classification.ipynb` and `R/main_figures/Fig5_TLS_classification_mfuzz.R` |
| Immunotherapy validation | `notebooks/main_figures/figure5-Immunotherapy.ipynb` |
| Metabolic analysis | `notebooks/main_figures/figure5-Metabolic Differences.ipynb` |
| Regulatory-network analysis | `notebooks/main_figures/figure5-pyscenic.ipynb` |
| Spatial multi-omic visualization | `notebooks/main_figures/figure6.ipynb` |
| COMMUNE workflow | `commune/` |

## Data availability

Human spatial and multi-omic data are deposited in the National Genomics Data Center under BioProjects `PRJCA068416` and `PRJCA068638`. Associated OMIX accessions are `OMIX018434`, `OMIX018461`, `OMIX018469` and `OMIX018470`. Access is controlled in accordance with the repository and applicable ethics requirements.

External immunotherapy datasets are available from the Gene Expression Omnibus under accessions `GSE243013` and `GSE169246`. Other public datasets used in the study are listed in the manuscript supplementary information.

## Project paths

The code does not contain machine-specific server paths. Set the `TLS_PROJECT_ROOT` environment variable to the directory containing the input data and analysis working directories.

Linux or macOS:

```bash
export TLS_PROJECT_ROOT=/path/to/project-data
```

Windows PowerShell:

```powershell
$env:TLS_PROJECT_ROOT = "D:\path\to\project-data"
```

If the variable is not set, scripts use the current working directory. Large data files and generated results are intentionally excluded from this repository.

## Reproducibility notes

Notebook outputs and execution counters have been cleared before public release. Run notebooks from the first cell so that `PROJECT_ROOT` is initialized. Software versions and detailed execution order should be added before the versioned archival release.

## Contact

For questions about the code or controlled-access data, please contact the corresponding authors listed in the manuscript.
