#!/bin/bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/rerender_convergence_%j.out
#SBATCH --error=scripts/figures/logs/rerender_convergence_%j.err

set -o pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

cd "${SCRIPT_DIR}"

echo "[$(date '+%H:%M:%S')] === fig5_convergence.R ==="
Rscript fig5_convergence.R || echo "WARN fig5_convergence.R FAILED"

echo "[$(date '+%H:%M:%S')] === fig5_convergence_v3.R ==="
Rscript fig5_convergence_v3.R || echo "WARN fig5_convergence_v3.R FAILED"

echo "[$(date '+%H:%M:%S')] === figS_hotspot_convergence.R ==="
Rscript figS_hotspot_convergence.R || echo "WARN figS_hotspot_convergence.R FAILED"

echo "[$(date '+%H:%M:%S')] === DONE ==="
