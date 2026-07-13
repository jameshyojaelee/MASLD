#!/bin/bash
# Protocol-clean LOO-dataset Hotspot re-runs for the 4 post-2026-05-22-remediation
# re-discovered cell types (mac/fib/chol/tcells). Mirrors the clean main run:
# always excludes the contaminated GSE136103 + Liver_Atlas (which also bumps the
# macrophage kNN to 50), and holds out each of the 5 clean datasets in turn.
# Hepatocytes are NOT re-run (coherent pre-remediation; legitimately use the
# Liver_Atlas hepatocyte fraction). 4 CT x 5 holdout = 20 tasks.
#
# Output: <ct>/loo/<holdout>/module_genes.tsv  (never touches the canonical full run)
# After all tasks: sbatch --dependency=afterok:<thisjob> 503b_loo_consolidate.sbatch
#SBATCH --job-name=Hotspot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=0-19
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/503_clean_%A_%a.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/503_clean_%A_%a.err

# Memory: the clean MAIN runs for these 4 non-global CTs peaked at 2.4-5.4 GB
# (only the global/1.2M-cell run needed ~170 GB), so 32 G is ample headroom and
# lets all 20 tasks run concurrently well under the nslab cap.
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

# Match the production/clean-run tuning.
export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=8

CELL_TYPES=(macrophages fibroblasts cholangiocytes tcells)
# 5 protocol-clean datasets (GSE136103 + Liver_Atlas excluded as contaminated).
DATASETS=(GSE174748 GSE185477 GSE189600 GSE202379 GSE244832)

CT_IDX=$(( SLURM_ARRAY_TASK_ID / 5 ))
DS_IDX=$(( SLURM_ARRAY_TASK_ID % 5 ))
CT=${CELL_TYPES[$CT_IDX]}
DS=${DATASETS[$DS_IDX]}

echo "[$(date)] CLEAN LOO task $SLURM_ARRAY_TASK_ID -> cell_type=$CT  hold-out=$DS  (exclude GSE136103+Liver_Atlas)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py \
  --cell-type "$CT" \
  --loo-dataset "$DS" \
  --exclude-datasets GSE136103 Liver_Atlas

echo "[$(date)] DONE task $SLURM_ARRAY_TASK_ID ($CT / hold-out $DS)"
