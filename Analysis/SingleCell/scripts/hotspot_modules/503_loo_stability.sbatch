#!/bin/bash
#SBATCH --job-name=hotspot_503
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-48%7
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/503_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/503_%A_%a.err

# LOO-dataset Hotspot stability — 49-task array (7 cell types × 7 datasets).
# %7 caps concurrency at 7 (one per cell type at a time); bigmem 512G matches
# the main 501 array config so hep-scale LOO tasks have headroom.
#
# Task index decomposition:
#   CT_IDX  = SLURM_ARRAY_TASK_ID / 7   (cell type)
#   DS_IDX  = SLURM_ARRAY_TASK_ID % 7   (dataset hold-out)
#
# After all 49 tasks finish, run the consolidator:
#   srun --partition=cpu --qos=interactive --mem=32G --cpus-per-task=4 --time=2:00:00 \
#     micromamba run -n rapids_singlecell python \
#     Analysis/SingleCell/scripts/hotspot_modules/503b_loo_consolidate.py

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

# Match the production array tuning (501_run_hotspot.sbatch).
export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=16

CELL_TYPES=(global hepatocytes macrophages fibroblasts endothelial_cells cholangiocytes tcells)
# Exact dataset strings from donor_metadata.tsv `dataset` column (enumerated 2026-05-13)
DATASETS=(GSE136103 GSE174748 GSE185477 GSE189600 GSE202379 GSE244832 Liver_Atlas)

CT_IDX=$(( SLURM_ARRAY_TASK_ID / 7 ))
DS_IDX=$(( SLURM_ARRAY_TASK_ID % 7 ))
CT=${CELL_TYPES[$CT_IDX]}
DS=${DATASETS[$DS_IDX]}

echo "[$(date)] LOO task $SLURM_ARRAY_TASK_ID -> cell_type=$CT  hold-out=$DS (max_genes=$HOTSPOT_MAX_GENES, jobs=$HOTSPOT_JOBS)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/503_loo_stability.py \
  --cell-type "$CT" \
  --dataset-holdout "$DS"

echo "[$(date)] DONE task $SLURM_ARRAY_TASK_ID"
