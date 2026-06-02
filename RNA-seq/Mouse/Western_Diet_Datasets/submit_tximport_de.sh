#!/bin/bash
#SBATCH --job-name=tximport_de
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/tximport_de_%j.out
#SBATCH --error=logs/tximport_de_%j.err

# =============================================================================
# Run tximport + limma-voom DE for a Western Diet dataset
# =============================================================================
# Usage:
#   sbatch submit_tximport_de.sh GSE220575
#   sbatch submit_tximport_de.sh GSE246088
#   sbatch submit_tximport_de.sh GSE305484
#
# Requires: Kallisto quant to be complete for the dataset
# =============================================================================

set -euo pipefail

DATASET_ID="${1:?ERROR: Supply DATASET_ID as first argument}"
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WDDIR="${BASE}/RNA-seq/Mouse/Western_Diet_Datasets"

echo "=== tximport + DE for ${DATASET_ID} ==="
echo "Date: $(date)"

# Activate rnaseq environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

# Verify Kallisto outputs exist
QUANT_DIR="${WDDIR}/${DATASET_ID}/quant/kallisto"
N_QUANT=$(find "${QUANT_DIR}" -name "abundance.tsv" 2>/dev/null | wc -l)
echo "Kallisto output directories with abundance.tsv: ${N_QUANT}"

if [[ ${N_QUANT} -eq 0 ]]; then
    echo "ERROR: No Kallisto outputs found in ${QUANT_DIR}"
    echo "  Run Kallisto quant first."
    exit 1
fi

# Run tximport + DE
Rscript "${WDDIR}/run_tximport_de.R" "${DATASET_ID}"

echo ""
echo "=== Completed: $(date) ==="
