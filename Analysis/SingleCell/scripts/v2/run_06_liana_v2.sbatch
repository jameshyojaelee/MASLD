#!/bin/bash
#SBATCH --job-name=06_liana_v2
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-34%4
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/06_liana_v2_%A_%a.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/06_liana_v2_%A_%a.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/per_donor_lr_v2

# Array task 1 also writes the lineage_cell_counts sidecar
EXTRA=""
if [ "${SLURM_ARRAY_TASK_ID:-0}" = "1" ]; then
  EXTRA="--write-counts"
fi

echo "[$(date)] starting 06_per_donor_liana_v2 array_task=${SLURM_ARRAY_TASK_ID} on $(hostname)"
micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/v2/06_per_donor_liana_v2.py \
    --array-task-id "${SLURM_ARRAY_TASK_ID}" \
    --array-batch-size 8 \
    ${EXTRA}
echo "[$(date)] DONE"
