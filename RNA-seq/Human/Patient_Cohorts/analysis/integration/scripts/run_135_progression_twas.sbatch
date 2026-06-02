#!/bin/bash
#SBATCH --job-name=135_prog_twas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# 135_progression_twas.R — Progression-specific TWAS re-analysis
#
# If existing S-PrediXcan TWAS results exist (ghodsian, chen), this only
# needs ~16 GB to load and cross-reference. If new TWAS runs are needed
# (FinnGen, BBJ), increase to bigmem/200G.
#
# Usage:
#   sbatch run_135_progression_twas.sbatch
#   # Or for a single GWAS:
#   GWAS_FILTER=ghodsian sbatch run_135_progression_twas.sbatch

set -euo pipefail

echo "=== Job ${SLURM_JOB_ID} on $(hostname) at $(date) ==="
echo "CPUs: ${SLURM_CPUS_PER_TASK}, Mem: ${SLURM_MEM_PER_NODE}MB"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

# Pass GWAS_FILTER if set
export GWAS_FILTER="${GWAS_FILTER:-ALL}"

Rscript "${SCRIPT_DIR}/135_progression_twas.R" 2>&1

echo "=== Job complete at $(date) ==="
