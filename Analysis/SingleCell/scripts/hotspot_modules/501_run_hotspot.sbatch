#!/bin/bash
#SBATCH --job-name=hotspot_501
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-6
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/501_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/501_%A_%a.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

# Tune: 2000 autocorr-significant genes (vs default 5000) gives ~6x speedup
# on local-corr step (O(genes^2 * cells)); 16 parallel jobs matches CPU count.
export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=16

CELL_TYPES=(global hepatocytes macrophages fibroblasts endothelial_cells cholangiocytes tcells)
CT=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}
echo "[$(date)] Array task $SLURM_ARRAY_TASK_ID -> $CT (max_genes=$HOTSPOT_MAX_GENES, jobs=$HOTSPOT_JOBS)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py \
  --cell-type "$CT"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/validate.py hotspot "$CT"

echo "[$(date)] Done $CT"
