#!/bin/bash
#SBATCH --job-name=dream_stage_contrasts
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/05c_dream_stage_contrasts_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/05c_dream_stage_contrasts_%j.err

set -euo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
cd "$SCRIPT_DIR"

echo "[start] $(date)"
echo "[host]  $(hostname)"
echo "[pwd]   $(pwd)"
echo "[cpus]  ${SLURM_CPUS_PER_TASK}"

micromamba run -n rnaseq Rscript 05c_dream_stage_contrasts.R

echo "[done] $(date)"
