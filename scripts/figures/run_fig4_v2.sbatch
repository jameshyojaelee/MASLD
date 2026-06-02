#!/bin/bash
#SBATCH --job-name=fig4_v2
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --output=logs/fig4_v2_%j.out
#SBATCH --error=logs/fig4_v2_%j.err

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$BASE"
mkdir -p logs

echo "[$(date)] Starting Figure 4 v2 build"
echo "Step 1/2: Pre-computation scan (evidence matrix)"
micromamba run -n rnaseq Rscript scripts/figures/fig4_module_evidence_scan.R

echo "[$(date)] Step 2/2: Main figure panels"
micromamba run -n rnaseq Rscript scripts/figures/fig4_validation_v2.R

echo "[$(date)] Done. Panels in figures/main/fig4_validation/panels/"
echo "Composite: figures/main/fig4_validation/fig4_validation.pdf"
