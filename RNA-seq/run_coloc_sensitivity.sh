#!/bin/bash
# run_coloc_sensitivity.sh
# Submits strand-ambiguity sensitivity (50b) and prior sensitivity (50c) as SLURM jobs.
# Power analysis (50d) runs interactively — instant computation.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== COLOC Sensitivity Analyses ==="
echo ""

# 50b: Strand-ambiguity sensitivity
JOB_50B=$(sbatch --parsable \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=32G \
  --time=4:00:00 \
  --job-name=coloc_strand_sens \
  --output="${SCRIPT_DIR}/logs/50b_strand_sensitivity_%j.out" \
  --error="${SCRIPT_DIR}/logs/50b_strand_sensitivity_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/50b_strand_ambiguity_sensitivity.R")
echo "  50b strand sensitivity: SLURM $JOB_50B"

# 50c: Prior sensitivity
JOB_50C=$(sbatch --parsable \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=16G \
  --time=2:00:00 \
  --job-name=coloc_prior_sens \
  --output="${SCRIPT_DIR}/logs/50c_prior_sensitivity_%j.out" \
  --error="${SCRIPT_DIR}/logs/50c_prior_sensitivity_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/50c_prior_sensitivity.R")
echo "  50c prior sensitivity: SLURM $JOB_50C"

echo ""
echo "  50d power analysis: run interactively (instant)"
echo "    micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/50d_power_analysis.R"
echo ""
echo "Jobs submitted. Monitor with: squeue -u \$USER"
