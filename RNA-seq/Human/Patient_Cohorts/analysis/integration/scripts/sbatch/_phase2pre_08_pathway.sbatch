#!/bin/bash
#SBATCH --job-name=p2pre_08
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p2pre_08_%j.out
#SBATCH --error=logs/p2pre_08_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)"
Rscript 08_pathway_analysis.R
echo "Finished: $(date)"
