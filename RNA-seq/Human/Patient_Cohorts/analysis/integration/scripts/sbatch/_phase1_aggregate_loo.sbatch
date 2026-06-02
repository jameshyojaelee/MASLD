#!/bin/bash
#SBATCH --job-name=p1_agg_loo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p1_agg_loo_%j.out
#SBATCH --error=logs/p1_agg_loo_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)"
Rscript aggregate_loo_cv.R
echo "Finished: $(date)"
