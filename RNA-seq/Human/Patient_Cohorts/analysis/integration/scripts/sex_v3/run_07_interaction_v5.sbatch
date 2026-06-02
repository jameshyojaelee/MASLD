#!/bin/bash
#SBATCH --job-name=v5_M07
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M07v5_%j.out
#SBATCH --error=logs/M07v5_%j.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
cd "${SCRIPT_DIR}"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
echo "=== Layer 1 interaction-test v5 start: $(date) ==="
Rscript 07_interaction_test_v5.R
echo "=== done: $(date) ==="
