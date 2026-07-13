# COMMUNE

COMMUNE identifies ecologically constrained multicellular communication modules
from annotated single-cell RNA-sequencing data. It combines sample-level
CellChat inference, NMF communication programs, and cell-type abundance
co-occurrence support.

## Workflow

1. Infer ligand-receptor communication independently for each sample with
   CellChat.
2. Build an interaction-by-sample probability matrix and decompose it with
   non-negative matrix factorization (NMF).
3. Estimate cell-type co-occurrence within each pathological group using
   Pearson correlation and Benjamini-Hochberg correction.
4. Retain NMF communication edges only when sender and receiver cell types have
   significant abundance-based co-occurrence support.

## System requirements

The tested environment is Linux x86_64 with Python 3.11.11 and R 4.4.3. No
non-standard hardware is required. The demo can run on a workstation with 8 GB
RAM; more memory and CPU cores are recommended for a full atlas.

Python dependencies are pinned in `requirements.txt`. Core R dependencies are:

- CellChat 2.2.0
- NMF 0.28
- Matrix 1.7.3
- dplyr 1.1.4
- tidyr 1.3.1
- tibble 3.2.1
- future 1.34.0

Check the active R environment with:

```bash
/usr/bin/Rscript scripts/check_dependencies.R
```

Create a Python environment with:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Installation normally takes 10-30 minutes when CellChat and its Bioconductor
dependencies are not already installed. CellChat installation instructions are
available from its upstream repository.

## Input data

The demo preparation script accepts an h5ad file whose `X` matrix contains raw
counts in CSR form. The following `obs` columns are required:

- `sample`: biological sample identifier
- `group`: pathological group
- `celltype_1`: broad immune/non-immune annotation
- `celltype_2`: intermediate cell-type annotation
- `celltype_3`: fine-grained cell-type annotation used by COMMUNE

The CellChat step consumes a genes-by-cells Matrix Market file plus matching
gene and cell metadata files. `prepare_demo_data.py` creates all of these from
the h5ad input and applies library-size normalization followed by `log1p` in the
CellChat script.

## Demo

Run the complete demo from the repository root:

```bash
bash run_demo.sh \
  /path/to/combined_adata_inner.h5ad \
  /path/to/commune_results
```

With no arguments, the demo reads:

```text
example_data/commune_demo.h5ad
```

and writes every intermediate result under:

```text
demo_results/
```

The demo deterministically selects 20 non-blood `Normal` samples and 20
non-blood `Tumor` samples, requiring both immune and non-immune compartments,
then retains at most 100 cells per sample. The demo uses rank 4, five NMF runs,
the five most abundant cell types per sample, and one CellChat bootstrap
iteration so that every stage can be exercised on a small dataset. The
co-occurrence thresholds remain Pearson `r > 0.3` and adjusted `P < 0.05`.
Reduced CellChat and NMF settings are for software validation, not manuscript
reproduction.

On the tested server, the complete demo took 212 seconds and used approximately
3.45 GB peak resident memory. Runtime on a typical desktop is expected to be
about 4-10 minutes, depending on disk speed and CPU performance.

Expected outputs are organized as:

```text
00_demo_data/              sampled h5ad, counts, genes, cells, and metadata
01_cooccurrence/           group-specific supported cell-type pairs
02_cellchat/               interaction matrices, long table, and sample QC
03_nmf/                    filtered matrices, NMF model, modules, assignments
04_constrained_modules/    constrained interactions, summaries, network plots
```

## Full analysis settings

For the manuscript analysis, use all eligible tissue samples and retain the
method parameters described in the paper:

- CellChat: `--min-interactions 500`, `--p-cutoff 0.05`
- Co-occurrence: `--correlation-cutoff 0.3`,
  `--adjusted-p-cutoff 0.05`
- NMF: `--min-samples 3`, `--top-features 10000`, `--rank 14`,
  `--nrun 200`; add `--survey-nrun 10` or more for rank assessment

Each script supports explicit input and output arguments, so full-data paths do
not need to be edited in source code. Run a script without arguments to see its
required options in the error message.

## Reproducibility notes

Random seeds are explicit in the demo data selection and NMF steps. Each R step
writes `sessionInfo.txt` beside its outputs. The long CellChat output preserves
sender, receiver, ligand, receptor, probability, P value, and sample fields so
that NMF feature identifiers can be audited without parsing ambiguous cell-type
names.

## License

COMMUNE is released under the MIT License. Third-party dependencies, including
CellChat, remain subject to their respective licenses.

## Citation

If you use COMMUNE, cite the accompanying manuscript. Repository
metadata are also provided in `CITATION.cff`.
