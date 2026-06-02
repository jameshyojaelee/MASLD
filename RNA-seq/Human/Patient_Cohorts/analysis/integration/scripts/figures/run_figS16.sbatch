#!/bin/bash
#SBATCH --job-name=figS16_plasma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS16_%j.out
#SBATCH --error=logs/figS16_%j.err

set -euo pipefail

echo "=== figS16: Plasma proteomics biological characterization ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/figures"
cd "$SCRIPT_DIR"

mkdir -p logs

Rscript figS16_plasma_proteomics.R

echo ""
echo "Completed: $(date)"
