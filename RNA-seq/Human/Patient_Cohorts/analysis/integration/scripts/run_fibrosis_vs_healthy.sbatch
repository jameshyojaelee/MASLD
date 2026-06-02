#!/bin/bash
#SBATCH --job-name=fib_vs_healthy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/fib_vs_healthy_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/fib_vs_healthy_%j.err
#
# Runs 05d_fibrosis_vs_healthy_dream.R to produce:
#   staging_classifier/fibrosis_vs_healthy_dream.csv
#
# Usage:
#   cd RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
#   sbatch run_fibrosis_vs_healthy.sbatch

set -euo pipefail

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd "${SCRIPT_DIR}"

echo "=== Fibrosis vs Healthy Dream ==="
echo "Host:  $(hostname)"
echo "Date:  $(date)"
echo "CPUs:  ${SLURM_CPUS_PER_TASK:-?}"
echo "PWD:   $(pwd)"
echo ""

micromamba run -n rnaseq Rscript "${SCRIPT_DIR}/05d_fibrosis_vs_healthy_dream.R"

echo "=== Done: $(date) ==="
