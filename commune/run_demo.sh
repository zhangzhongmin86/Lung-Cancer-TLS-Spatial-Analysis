#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_H5AD="${1:-${ROOT_DIR}/example_data/commune_demo.h5ad}"
RESULTS_ROOT="${2:-${ROOT_DIR}/demo_results}"
R_SCRIPT="${R_SCRIPT:-/usr/bin/Rscript}"

DEMO_DIR="${RESULTS_ROOT}/00_demo_data"
COOCCURRENCE_DIR="${RESULTS_ROOT}/01_cooccurrence"
CELLCHAT_DIR="${RESULTS_ROOT}/02_cellchat"
NMF_DIR="${RESULTS_ROOT}/03_nmf"
MODULE_DIR="${RESULTS_ROOT}/04_constrained_modules"

mkdir -p "${RESULTS_ROOT}"

python "${ROOT_DIR}/scripts/prepare_demo_data.py" \
  --input "${INPUT_H5AD}" \
  --output-dir "${DEMO_DIR}" \
  --groups Normal Tumor \
  --samples-per-group 20 \
  --cells-per-sample 100 \
  --max-celltypes-per-sample 12 \
  --seed 20250713

"${R_SCRIPT}" "${ROOT_DIR}/Abundance-based cell-type co-occurrence estimation.R" \
  --metadata "${DEMO_DIR}/metadata.csv" \
  --output-dir "${COOCCURRENCE_DIR}" \
  --groups Normal,Tumor \
  --min-samples 3 \
  --correlation-cutoff 0.3 \
  --adjusted-p-cutoff 0.05

"${R_SCRIPT}" "${ROOT_DIR}/Ligand–receptor interaction inference.R" \
  --counts "${DEMO_DIR}/counts.mtx.gz" \
  --genes "${DEMO_DIR}/genes.tsv" \
  --metadata "${DEMO_DIR}/metadata.csv" \
  --output-dir "${CELLCHAT_DIR}" \
  --groups Normal,Tumor \
  --min-cells 3 \
  --max-celltypes 5 \
  --min-interactions 1 \
  --p-cutoff 0.05 \
  --nboot 1 \
  --seed 1

"${R_SCRIPT}" "${ROOT_DIR}/NMF-based communication program extraction.R" \
  --matrix "${CELLCHAT_DIR}/cellchat_probability_interaction_by_sample.csv" \
  --interactions "${CELLCHAT_DIR}/cellchat_interactions_long.csv" \
  --output-dir "${NMF_DIR}" \
  --min-samples 2 \
  --top-features 500 \
  --rank 4 \
  --nrun 5 \
  --survey-nrun 0 \
  --seed 100

python "${ROOT_DIR}/scripts/construct_modules.py" \
  --modules "${NMF_DIR}/nmf_module_interactions.csv" \
  --cooccurrence-dir "${COOCCURRENCE_DIR}" \
  --metadata "${DEMO_DIR}/metadata.csv" \
  --output-dir "${MODULE_DIR}" \
  --top-k 20

echo "COMMUNE demo completed: ${RESULTS_ROOT}"
