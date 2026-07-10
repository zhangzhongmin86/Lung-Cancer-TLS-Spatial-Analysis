# Lung-Cancer-TLS-Spatial-Analysis

Analysis code for the manuscript **“A tumor-proximal EPAS1-high lipid-stress niche diverts tertiary lymphoid structure maturation in lung adenocarcinoma.”**

## Repository contents

- `figure*.R`: R scripts used for main and supplementary figures.
- `figure*.ipynb`: Python notebooks used for spatial, metabolic, immunotherapy and regulatory-network analyses.
- `commune/`: scripts and notebooks implementing the COMMUNE workflow.

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
