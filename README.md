# Lung-Cancer-TLS-Spatial-Analysis

Analysis code for the accompanying manuscript.

## Repository contents

- `R/main_figures/`: R script for Figure 1.
- `R/supplementary_figures/`: R scripts for supplementary figures.
- `notebooks/`: Python notebooks.
- `commune/`: COMMUNE workflow implementation.

## Main analysis index

| Analysis | Code |
| --- | --- |
| Atlas construction | `R/main_figures/figure1.R` |
| COMMUNE workflow | `commune/` |
| TLS identification | `notebooks/main_figures/figure3-Identifying_TLS.ipynb` |
| Regional grouping | `notebooks/main_figures/figure3-Regional_grouping.ipynb` |
| Cell-subtype subclustering (S1-S8) | `R/supplementary_figures/figureS1.R` - `figureS8.R` |
| NMF cell composition (S9-S10) | `R/supplementary_figures/figureS9_S10.R` |
| EPAS1 spatial features (S29) | `notebooks/supplementary_figures/figureS29.ipynb` |
| Reaction-step validation (S30) | `notebooks/supplementary_figures/figureS30.ipynb` |

## Data availability

Human spatial and multi-omic data are deposited in the National Genomics Data Center under BioProjects `PRJCA068416` and `PRJCA068638`. Associated OMIX accessions are `OMIX018434`, `OMIX018461`, `OMIX018469` and `OMIX018470`.

External immunotherapy datasets: `GSE243013`, `GSE169246`.

## Project paths

Set `TLS_PROJECT_ROOT` environment variable to the directory containing input data.

## License

See `commune/LICENSE`.
