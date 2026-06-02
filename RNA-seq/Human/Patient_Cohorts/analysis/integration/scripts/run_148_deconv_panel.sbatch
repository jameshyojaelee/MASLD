#!/bin/bash
#SBATCH --job-name=s148_deconv_panel
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/148_deconv_gene_panel_%j.out
#SBATCH --error=logs/148_deconv_gene_panel_%j.err

set -euo pipefail

echo "=== Job $SLURM_JOB_ID started at $(date) ==="
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
mkdir -p "${SCRIPT_DIR}/../logs"

Rscript "${SCRIPT_DIR}/148_deconv_gene_panel.R"

echo "=== Job completed at $(date) ==="
