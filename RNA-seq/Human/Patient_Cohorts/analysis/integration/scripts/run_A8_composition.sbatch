#!/bin/bash
#SBATCH --job-name=A8_compadj
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/A8_compadj_%j.out
#SBATCH --error=logs/A8_compadj_%j.err

# A8: Composition-adjusted sex-interaction dream refit
# Three model variants on the 5-cohort canonical:
#   M0 = canonical (baseline reproduction)
#   M1 = +hep_frac
#   M2 = +hep_frac + mac_frac + endo_frac + chol_frac
# Headline output: composition_independent_genes.csv

set -eo pipefail
mkdir -p logs

echo "=== A8 composition-adjusted sex-interaction refit ==="
echo "Job ID:   ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Started:  $(date)"
echo "CPUs:     ${SLURM_CPUS_PER_TASK}"
echo "Mem:      ${SLURM_MEM_PER_NODE} MB"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42

Rscript 14.4f_composition_adjusted_T13.R

echo "Done: $(date)"
