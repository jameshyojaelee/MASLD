#!/bin/bash
#SBATCH --job-name=Hotspot
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-5
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/501_clean_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/501_clean_%A_%a.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=16

# 6 non-hepatocyte cell types; hepatocytes skipped (NPC datasets contribute essentially 0 hep cells, no contamination)
CELL_TYPES=(global macrophages fibroblasts endothelial_cells cholangiocytes tcells)
CT=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}
echo "[$(date)] Array task $SLURM_ARRAY_TASK_ID -> $CT (clean run, excluding GSE136103+Liver_Atlas)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py \
  --cell-type "$CT" \
  --exclude-datasets GSE136103 Liver_Atlas

echo "[$(date)] Done $CT"
