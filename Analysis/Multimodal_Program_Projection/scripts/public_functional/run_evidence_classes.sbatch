#!/bin/bash
#SBATCH --job-name=pubclass
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/class_%A_%a.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/class_%A_%a.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/public_functional"
DATASETS=(GSE200418 GSE253380 GSE106737)
cd "${PROJECT_ROOT}"
micromamba run -n rnaseq Rscript "${SCRIPT_ROOT}/07_run_evidence_classes.R" "${DATASETS[${SLURM_ARRAY_TASK_ID}]}"
