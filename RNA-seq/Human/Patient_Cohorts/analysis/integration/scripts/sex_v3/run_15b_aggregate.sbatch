#!/bin/bash
#SBATCH --job-name=v6_P6agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P6agg_%j.out
#SBATCH --error=logs/P6agg_%j.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== P6 aggregator ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 15b_cochranQ_aggregate.R
echo "=== done: $(date) ==="
