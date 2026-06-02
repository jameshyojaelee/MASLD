#!/bin/bash
#SBATCH --job-name=loo_sensitivity_reexec
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/loo_sensitivity_reexec_%j.log

# Re-run LOO-CV sensitivity with new DEG thresholds (padj<0.05, |logFC|>0.3)

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${BASE}"

echo "=== LOO sensitivity re-execution started: $(date) ==="

cd RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
Rscript 05c_gse135251_loo_sensitivity.R 2>&1

echo "=== Completed: $(date) ==="
