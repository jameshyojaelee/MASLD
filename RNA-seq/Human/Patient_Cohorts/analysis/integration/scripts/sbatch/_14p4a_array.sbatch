#!/bin/bash
#SBATCH --job-name=14p4a_arr
#SBATCH --array=1-50
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/14p4a_arr_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/14p4a_arr_%A_%a.err

# 14.4a_sex_permutation_iter array. Each task = 1 permutation.

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export PERM_ITER=${SLURM_ARRAY_TASK_ID}

echo "Started: $(date)  PERM_ITER=${PERM_ITER}  CPUs=${SLURM_CPUS_PER_TASK}"
Rscript 14.4a_sex_permutation_iter.R
echo "Finished: $(date)"
