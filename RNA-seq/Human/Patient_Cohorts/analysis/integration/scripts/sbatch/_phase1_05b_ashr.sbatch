#!/bin/bash
#SBATCH --job-name=p1_05b_ashr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p1_05b_ashr_%j.out
#SBATCH --error=logs/p1_05b_ashr_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)"
Rscript 05b_ashr_shrinkage.R
echo "Finished: $(date)"
