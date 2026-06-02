#!/usr/bin/env bash
#SBATCH --job-name=v1v2_concord
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2_phase05/logs/v1v2_concord_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2_phase05/logs/v1v2_concord_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/Analysis/SingleCell/scripts/v2/08_v1_vs_v2_concordance.R"

mkdir -p "${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2_phase05/logs"
mkdir -p "${PROJECT_ROOT}/figures/supplementary/stage_ccc_v2"

echo "[start] $(date)"
echo "[host]  $(hostname)"
echo "[script] ${SCRIPT}"

micromamba run -n rnaseq Rscript "${SCRIPT}"

echo "[done] $(date)"
