#!/bin/bash
#SBATCH --job-name=sex_T4_boot
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=1-40%6
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap/logs/T4_boot_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap/logs/T4_boot_%A_%a.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42
export BOOT_ITER=${SLURM_ARRAY_TASK_ID}
export N_BOOT=5

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript ${PROJ}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.4e_sex_bootstrap_T4.R
