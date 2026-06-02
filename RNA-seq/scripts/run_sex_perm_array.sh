#!/bin/bash
#SBATCH --array=1-100%8
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --qos=nslab
#SBATCH --job-name=permutation
#SBATCH --output=RNA-seq/scripts/logs/sex_perm_%A_%a.out
#SBATCH --error=RNA-seq/scripts/logs/sex_perm_%A_%a.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/scripts/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

export TASK_ID=$SLURM_ARRAY_TASK_ID
export N_PER_TASK=100

Rscript RNA-seq/scripts/14_4_sex_perm_array_task.R
