# Example data

The bundled demonstration subset contains 40 samples (20 `Normal` and 20
`Tumor`) and 4,000 cells in `commune_demo.h5ad`. It was generated reproducibly
with:

```bash
python scripts/prepare_demo_data.py \
  --input /path/to/combined_adata_inner.h5ad \
  --output-dir example_data/generated \
  --groups Normal Tumor \
  --samples-per-group 20 \
  --cells-per-sample 100 \
  --max-celltypes-per-sample 12 \
  --seed 20250713
```

The bundled h5ad file is the default input used by `run_demo.sh`. Intermediate
Matrix Market files and all downstream outputs are generated under
`demo_results/` and are not stored in the repository.
