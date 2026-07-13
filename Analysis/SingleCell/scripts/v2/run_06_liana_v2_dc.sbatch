#!/bin/bash
#SBATCH --job-name=liana
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-16
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/06_liana_v2_dc_%A_%a.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/06_liana_v2_dc_%A_%a.err

# Donor-collapsed LIANA (pseudoreplication fix). 121 true donors / batch 8 = 16
# tasks. Pools cells per biological donor and runs LIANA once per donor.
# Outputs -> per_donor_lr_v2_dc/ (run-level per_donor_lr_v2/ left intact).
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/per_donor_lr_v2_dc

# Array task 1 also writes the donor x lineage cell-count sidecar
EXTRA=""
if [ "${SLURM_ARRAY_TASK_ID:-0}" = "1" ]; then
  EXTRA="--write-counts"
fi

echo "[$(date)] starting 06_per_donor_liana_v2 (DONOR-COLLAPSE) array_task=${SLURM_ARRAY_TASK_ID} on $(hostname)"
micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/v2/06_per_donor_liana_v2.py \
    --array-task-id "${SLURM_ARRAY_TASK_ID}" \
    --array-batch-size 8 \
    ${EXTRA}
echo "[$(date)] DONE task ${SLURM_ARRAY_TASK_ID}"
