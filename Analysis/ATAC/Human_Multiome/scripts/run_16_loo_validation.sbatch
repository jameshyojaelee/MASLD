#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/16_loo_%j.out
#SBATCH --error=logs/16_loo_%j.err

set -eo pipefail

# -----------------------------------------------------------------------------
# B4-Fork-1: internal LOO + jackknife validation of GSE244832 scATAC
#
# Replacement for an external scATAC replication cohort (GSE296875 was confirmed
# unsuitable — no MASLD disease staging). Quantifies stability of:
#   - the 24 hepatocyte disease regulons (LOO over 18 donors + 10/25/50%
#     cell-subset jackknife)
#   - chromVAR top-20 TFs per cell type (LOO over 18 donors)
#
# bigmem + interactive QOS bypasses the 7TB nslab memory cap. Only one such job
# at a time (interactive QOS limits to 4 concurrent jobs). HDF5_USE_FILE_LOCKING
# is unset per B2 lesson learned (parallel array tasks share the h5 file).
# -----------------------------------------------------------------------------

echo "=== B4-Fork-1 LOO/jackknife validation ==="
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      $(hostname)"
echo "CPUs:      ${SLURM_CPUS_PER_TASK}"
echo "Memory:    ${SLURM_MEM_PER_NODE}MB"
echo "Started:   $(date)"
echo ""

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export HDF5_USE_FILE_LOCKING=FALSE

WORK_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome"
cd "${WORK_DIR}"
mkdir -p logs results/loo_validation

python scripts/16_loo_donor_validation.py

echo ""
echo "=== Verifying outputs ==="
OUT="${WORK_DIR}/results/loo_validation"
for f in loo_regulon_per_donor.csv loo_regulon_summary.csv \
         jackknife_regulon.csv jackknife_regulon_summary.csv \
         loo_chromvar_per_celltype.csv loo_chromvar_summary.csv \
         loo_outlier_flag.csv loo_validation_report.txt; do
    full="${OUT}/${f}"
    if [[ -f "${full}" ]]; then
        lines=$(wc -l < "${full}")
        echo "  OK: ${f} (${lines} lines)"
    else
        echo "  MISSING: ${f}"
    fi
done

echo ""
echo "Finished: $(date)"
