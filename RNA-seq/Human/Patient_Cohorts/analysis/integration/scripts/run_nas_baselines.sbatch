#!/bin/bash
#SBATCH --job-name=nas_baselines
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/nas_baselines_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/nas_baselines_%j.err
#
# Runs 05c_nas_vs_baseline_dream.R to produce:
#   nas_vs_nas0_dream.csv, nas_vs_healthy_dream.csv
#
# Usage:
#   cd RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
#   sbatch run_nas_baselines.sbatch

set -euo pipefail

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd "${SCRIPT_DIR}"

echo "=== NAS vs Baseline Dream ==="
echo "Host:  $(hostname)"
echo "Date:  $(date)"
echo "CPUs:  ${SLURM_CPUS_PER_TASK:-?}"
echo "PWD:   $(pwd)"
echo ""

micromamba run -n rnaseq Rscript "${SCRIPT_DIR}/05c_nas_vs_baseline_dream.R"

echo "=== Done: $(date) ==="
