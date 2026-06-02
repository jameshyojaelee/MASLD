#!/bin/bash
#SBATCH --job-name=ssgsea_fix
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/logs/ssgsea_fix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/logs/ssgsea_fix_%j.err

set -e

PATHWAY_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Pathway"
SCRIPTS_DIR="${PATHWAY_DIR}/scripts"
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PREFIX="${PATHWAY_DIR}/.mamba/pathway_analysis"

echo "=== Running Fixed ssGSEA Pipeline ==="
echo "Timestamp: $(date)"

# Step 1: Run refactored ssGSEA (with fixes)
echo ""
echo "[1/4] Running ssGSEA (no subsetting)..."
${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/run_refactored_ssgsea.R" 2>&1

# Step 2: Run advanced plots
echo ""
echo "[2/4] Generating advanced plots..."
${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/generate_advanced_plots.R" 2>&1

# Step 3: Run species contrast plots
echo ""
echo "[3/4] Generating species contrast plots..."
${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/generate_species_contrast_plots.R" 2>&1

# Step 4: Run creative plots
echo ""
echo "[4/4] Generating creative plots..."
${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/generate_creative_plots.R" 2>&1

echo ""
echo "=== Pipeline Complete ==="
echo "Output plots in: ${PATHWAY_DIR}/plots/refactored/"
echo "                 ${PATHWAY_DIR}/plots/creative/"
