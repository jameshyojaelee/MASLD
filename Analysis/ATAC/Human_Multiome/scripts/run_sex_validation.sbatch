#!/bin/bash
#SBATCH --job-name=chromVAR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/sex_validation_%j.out
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/sex_validation_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
eval "$(micromamba shell hook -s bash)"

ATAC_DIR="Analysis/ATAC/Human_Multiome"
RESULTS="${ATAC_DIR}/results/sex_validation"
mkdir -p "${RESULTS}" "${ATAC_DIR}/logs"

echo "=== Step 1: Sex epigenomic validation (Python, snapatac2 env) ==="
echo "Start: $(date)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1

python "${ATAC_DIR}/scripts/13_sex_epigenomic_validation.py" \
    --output-dir "${RESULTS}"

echo ""
echo "=== Step 2: Generate figures (R, rnaseq env) ==="
echo "Start: $(date)"
micromamba activate rnaseq

Rscript scripts/figures/figS_sex_epigenomic_validation.R 2>&1

echo ""
echo "=== Done: $(date) ==="
echo "Results:"
ls -la "${RESULTS}/"
echo "Figures:"
ls -la figures/supplementary/epigenomic/panel_sex_validation_*.pdf 2>/dev/null || echo "  (no figures yet)"
