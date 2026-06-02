#!/bin/bash
#SBATCH --job-name=15c_arr
#SBATCH --array=1-60
#SBATCH --partition=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/15c_arr_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/15c_arr_%A_%a.err

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)  ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}  CPUs=${SLURM_CPUS_PER_TASK}"
Rscript 15c_one_iter.R
echo "Finished: $(date)"
