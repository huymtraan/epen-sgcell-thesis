#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=epen-copykat
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/copykat_%j.out
#SBATCH --error=logs/copykat_%j.err

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${EPEN_OUTPUT_DIR:-${REPO_DIR}/outputs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
mkdir -p "${REPO_DIR}/logs" "${OUTPUT_DIR}/copykat_per_sample"
"${PYTHON_BIN}" "${SCRIPT_DIR}/06_copykat.py" \
  --input "${OUTPUT_DIR}/04_clustered.h5ad" \
  --output-dir "${OUTPUT_DIR}/copykat_per_sample" \
  --threads "${SLURM_CPUS_PER_TASK:-8}"
