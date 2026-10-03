#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=epen-processing
#SBATCH --cpus-per-task=36
#SBATCH --mem=120G
#SBATCH --output=logs/processing_%j.out
#SBATCH --error=logs/processing_%j.err

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${SCRIPT_DIR}/../logs"
bash "${SCRIPT_DIR}/run_processing.sh"
