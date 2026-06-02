#!/bin/bash
#SBATCH --job-name=v6_P5agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P5agg_%j.out
#SBATCH --error=logs/P5agg_%j.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== P5 aggregator ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 14b_calibration_summary.R
echo "=== done: $(date) ==="
