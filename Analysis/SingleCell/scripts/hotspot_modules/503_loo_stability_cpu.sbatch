#!/bin/bash
#SBATCH --job-name=hotspot_503cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=4-00:00:00
#SBATCH --array=3-48
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/503cpu_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/503cpu_%A_%a.err

# LOO-dataset Hotspot stability — re-submission for pending tasks 3-48
# (tasks 0,1,2 = global LOO on GSE136103/GSE174748/GSE185477 are still running
# on bigmem under job 16053346). No %N concurrency cap — let them run wild;
# scheduler will fit as many as cluster has free.
#
# cpu partition: 32 cpus/200G is ample (empirical MaxRSS for chol smoke was
# 3 GB; global with 1.2M cells extrapolates to ~100 GB, well under 200 GB).
# Bumping HOTSPOT_JOBS=32 to match the new cpu count.

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=32

CELL_TYPES=(global hepatocytes macrophages fibroblasts endothelial_cells cholangiocytes tcells)
DATASETS=(GSE136103 GSE174748 GSE185477 GSE189600 GSE202379 GSE244832 Liver_Atlas)

CT_IDX=$(( SLURM_ARRAY_TASK_ID / 7 ))
DS_IDX=$(( SLURM_ARRAY_TASK_ID % 7 ))
CT=${CELL_TYPES[$CT_IDX]}
DS=${DATASETS[$DS_IDX]}

echo "[$(date)] LOO task $SLURM_ARRAY_TASK_ID -> cell_type=$CT  hold-out=$DS (max_genes=$HOTSPOT_MAX_GENES, jobs=$HOTSPOT_JOBS, partition=cpu)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/503_loo_stability.py \
  --cell-type "$CT" \
  --dataset-holdout "$DS"

echo "[$(date)] DONE task $SLURM_ARRAY_TASK_ID"
