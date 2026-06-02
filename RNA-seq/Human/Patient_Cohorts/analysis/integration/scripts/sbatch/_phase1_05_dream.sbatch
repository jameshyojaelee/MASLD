#!/bin/bash
#SBATCH --job-name=p1_05_dream
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p1_05_dream_%j.out
#SBATCH --error=logs/p1_05_dream_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Host: $(hostname)"
echo "Started: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript 05_dream_mega_analysis.R

echo "Finished: $(date)"
