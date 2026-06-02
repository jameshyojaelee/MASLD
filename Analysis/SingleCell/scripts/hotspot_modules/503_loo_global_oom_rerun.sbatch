#!/bin/bash
#SBATCH --job-name=hotspot_503rerun
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=4-00:00:00
#SBATCH --array=3,4
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/503rerun_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/503rerun_%A_%a.err

# Rerun of cpu LOO indices 3 and 4 that hit OUT_OF_MEMORY at 200G:
#   idx=3 -> CT=global, DS=GSE189600
#   idx=4 -> CT=global, DS=GSE202379
# Global Hotspot LOO needs ~100-150 GB working set on 1.0M-cell subsets.

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=16

CELL_TYPES=(global hepatocytes macrophages fibroblasts endothelial_cells cholangiocytes tcells)
DATASETS=(GSE136103 GSE174748 GSE185477 GSE189600 GSE202379 GSE244832 Liver_Atlas)

CT_IDX=$(( SLURM_ARRAY_TASK_ID / 7 ))
DS_IDX=$(( SLURM_ARRAY_TASK_ID % 7 ))
CT=${CELL_TYPES[$CT_IDX]}
DS=${DATASETS[$DS_IDX]}

echo "[$(date)] LOO rerun task $SLURM_ARRAY_TASK_ID -> cell_type=$CT  hold-out=$DS (bigmem 512G)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/503_loo_stability.py \
  --cell-type "$CT" \
  --dataset-holdout "$DS"

echo "[$(date)] DONE task $SLURM_ARRAY_TASK_ID"
