#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${EPEN_OUTPUT_DIR:-${REPO_DIR}/outputs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
: "${EPEN_INPUT_H5AD:?Set EPEN_INPUT_H5AD to the concatenated raw-count H5AD}"

mkdir -p "${OUTPUT_DIR}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/01_qc.py" \
  --input "${EPEN_INPUT_H5AD}" --output "${OUTPUT_DIR}/01_qc.h5ad" \
  --nmads 3 --min-genes 100 --min-cells 3 --threads "${SLURM_CPUS_PER_TASK:-8}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/02_normalize_pca.py" \
  --input "${OUTPUT_DIR}/01_qc.h5ad" --output "${OUTPUT_DIR}/02_norm_pca.h5ad" \
  --n-top-genes 2000 --n-pcs 50 --threads "${SLURM_CPUS_PER_TASK:-8}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/03_harmony.py" \
  --input "${OUTPUT_DIR}/02_norm_pca.h5ad" --output "${OUTPUT_DIR}/03_harmony.h5ad" \
  --n-neighbors 15 --threads "${SLURM_CPUS_PER_TASK:-8}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/04_clustering.py" \
  --input "${OUTPUT_DIR}/03_harmony.h5ad" --output "${OUTPUT_DIR}/04_clustered.h5ad" \
  --threads "${SLURM_CPUS_PER_TASK:-8}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/05_differential_expression.py" \
  --input "${OUTPUT_DIR}/04_clustered.h5ad" --output "${OUTPUT_DIR}/04_clustered.h5ad"
