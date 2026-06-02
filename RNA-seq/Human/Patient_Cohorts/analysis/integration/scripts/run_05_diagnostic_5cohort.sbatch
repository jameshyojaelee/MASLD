#!/bin/bash
#SBATCH --job-name=dream_diag_5co
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/dream_diag_5co_%j.out
#SBATCH --error=logs/dream_diag_5co_%j.err

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Host: $(hostname)"
echo "Started: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript 05_diagnostic_5cohort_re_vs_fe.R

echo "Finished: $(date)"
