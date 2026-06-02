#!/bin/bash
#SBATCH --job-name=sex_T1_pmsub
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-50%10
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_power_matched/logs/T1_pm_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_power_matched/logs/T1_pm_%A_%a.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42
export SUBSAMPLE_ITER=${SLURM_ARRAY_TASK_ID}
export N_PER_CELL_M=52

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript ${PROJ}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.4d_sex_power_matched_T1.R
