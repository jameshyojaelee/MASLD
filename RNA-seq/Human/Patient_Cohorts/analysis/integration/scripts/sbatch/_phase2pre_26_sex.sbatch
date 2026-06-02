#!/bin/bash
#SBATCH --job-name=p2pre_26
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p2pre_26_%j.out
#SBATCH --error=logs/p2pre_26_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)"
Rscript 26_sex_stratified_analysis.R
echo "Finished: $(date)"
