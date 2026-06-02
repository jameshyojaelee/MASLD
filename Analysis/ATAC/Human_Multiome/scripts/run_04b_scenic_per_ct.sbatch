#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=72:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/04b_scenic_%A_%a.out
#SBATCH --error=logs/04b_scenic_%A_%a.err

set -eo pipefail

# -----------------------------------------------------------------------------
# B2 ATAC improvement plan: per-cell-type SCENIC+ regulons
#
# Runs `04b_scenic_plus_grn_per_celltype.py` once per array task on four
# non-hepatocyte populations so that C4 Panel J (multi-CT disease regulon
# heatmap) has analogous inputs to the existing hepatocyte set.
#
# Array tasks: Macrophage+Kupffer, Stellate, Endothelial+LSEC, Cholangiocyte
# bigmem + interactive QOS bypasses the 7TB nslab memory cap (interactive QOS
# allows 4 concurrent jobs, matching the array size). 72h wall is generous;
# Script 04b finished hepatocytes in ~4 min on 15K cells.
# -----------------------------------------------------------------------------

CELLTYPES=(
  "Macrophage,Kupffer_Cell:Macrophage"
  "Stellate_Cell:Stellate"
  "Endothelial,LSEC:Endothelial"
  "Cholangiocyte:Cholangiocyte"
)
PAIR=${CELLTYPES[$SLURM_ARRAY_TASK_ID]}
LABELS=${PAIR%:*}
NAME=${PAIR#*:}

echo "=== B2 per-CT SCENIC+ (array task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Job ID:    ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node:      $(hostname)"
echo "Cell type: ${NAME}"
echo "Labels:    ${LABELS}"
echo "CPUs:      ${SLURM_CPUS_PER_TASK}"
echo "Memory:    ${SLURM_MEM_PER_NODE}MB"
echo "Started:   $(date)"
echo ""

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
# Allow concurrent read of the same H5 dataset by parallel array tasks
# (fix for tasks 0/2 RuntimeError: H5Fopen() unable to lock file, errno=11)
export HDF5_USE_FILE_LOCKING=FALSE

WORK_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome"
cd "${WORK_DIR}"
mkdir -p logs scenic_plus

python scripts/04b_scenic_plus_grn_per_celltype.py \
    --target-labels "${LABELS}" \
    --celltype-name "${NAME}"

echo ""
echo "=== Verifying outputs for ${NAME} ==="
TAG=$(echo "${NAME}" | tr '[:upper:]' '[:lower:]')
for f in "${TAG}_regulons.csv" "${TAG}_enhancer_gene_links.csv" \
         "${TAG}_regulon_activity_scores.csv" "${TAG}_disease_regulons.csv"; do
    full="${WORK_DIR}/scenic_plus/${f}"
    if [[ -f "${full}" ]]; then
        lines=$(wc -l < "${full}")
        echo "  OK: ${f} (${lines} lines)"
    else
        echo "  MISSING: ${f}"
    fi
done

echo ""
echo "Finished: $(date)"
