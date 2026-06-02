#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=1-22
#SBATCH --job-name=susie_coloc
#SBATCH --output=logs/susie_coloc_%x_%A_%a.out
#SBATCH --error=logs/susie_coloc_%x_%A_%a.err

# SuSiE-COLOC array job — one chromosome per task
# Supports N_WORKERS>1 for within-node parallelization.
# Usage: GWAS_NAME=UKBB_ALT sbatch src/06_susie_coloc.sh
#        GWAS_NAME=UKBB_ALT N_WORKERS=8 sbatch --cpus-per-task=8 --mem=96G src/06_susie_coloc.sh

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

export CHR_FILTER="${SLURM_ARRAY_TASK_ID}"
N_WORKERS="${N_WORKERS:-1}"

echo "=== SuSiE-COLOC: ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "Workers: ${N_WORKERS} | CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Start: $(date)"

if [ "${N_WORKERS}" -le 1 ]; then
  # Single worker — run directly
  Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}"
else
  # Multi-worker: launch N parallel R processes with stride pattern
  PIDS=()
  for w in $(seq 0 $((N_WORKERS - 1))); do
    WORKER_ID=${w} N_WORKERS=${N_WORKERS} \
      Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}" &
    PIDS+=($!)
  done
  echo "Launched ${N_WORKERS} workers. PIDs: ${PIDS[*]}"

  # Wait for all workers
  FAIL=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=$((FAIL + 1))
  done
  if [ $FAIL -eq $N_WORKERS ]; then
    echo "ERROR: ALL ${N_WORKERS} workers failed — aborting"
    exit 1
  elif [ $FAIL -gt 0 ]; then
    echo "WARNING: ${FAIL}/${N_WORKERS} workers failed — merging partial results"
  fi

  # Merge worker outputs into final file
  OUT_DIR="results/susie_coloc/${GWAS_NAME}"
  OUT_FILE="${OUT_DIR}/susie_coloc_chr${CHR_FILTER}.csv"
  WORKER_FILES=$(ls ${OUT_DIR}/worker_chr${CHR_FILTER}_w*.csv 2>/dev/null || true)

  if [ -n "$WORKER_FILES" ]; then
    echo "Merging worker outputs..."
    Rscript -e "
      library(data.table)
      files <- commandArgs(trailingOnly = TRUE)
      dt_list <- lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL))
      dt_list <- dt_list[!sapply(dt_list, is.null)]
      if (length(dt_list) > 0) {
        merged <- rbindlist(dt_list, fill = TRUE)
        merged <- merged[!duplicated(ensembl)]
        fwrite(merged, '${OUT_FILE}')
        cat('Merged', nrow(merged), 'genes into final output\n')
      }
    " ${WORKER_FILES}
    # Cleanup worker temp files
    rm -f ${OUT_DIR}/worker_chr${CHR_FILTER}_w*.csv 2>/dev/null || true
    rm -f ${OUT_DIR}/checkpoint_chr${CHR_FILTER}_w*.csv 2>/dev/null || true
  fi
fi

# Cleanup shared checkpoint
rm -f "results/susie_coloc/${GWAS_NAME}/checkpoint_chr${CHR_FILTER}.csv" 2>/dev/null || true

echo "End: $(date)"
