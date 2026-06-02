#!/bin/bash
#SBATCH --array=1-10000%50
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --job-name=CPSS
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/cpss_%A_%a.out
#SBATCH --error=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/cpss_%A_%a.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

export CPSS_PAIR=$SLURM_ARRAY_TASK_ID
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/dream_cpss_iter.R
