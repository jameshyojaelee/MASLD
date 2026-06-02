#!/usr/bin/env bash
#SBATCH --job-name=343p_bulk_gsea_LR
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343p_bulk_gsea_LR_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343p_bulk_gsea_LR_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/Analysis/SingleCell/scripts/343p_bulk_gsea_LR_signatures.R"

mkdir -p "${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

micromamba run -n rnaseq Rscript "${SCRIPT}"
