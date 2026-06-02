#!/bin/bash -l
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=210G
#SBATCH --time=90:00:00  # capped at 90h per CLAUDE.md absolute-max policy (T2.3 2026-04-22); per-gene stride + resume handles re-runs
#SBATCH --job-name=susie_coloc
#SBATCH --output=logs/susie_coloc_bigmem_%A_%a.out
#SBATCH --error=logs/susie_coloc_bigmem_%A_%a.err

# Bigmem SuSiE-COLOC: runs N_WORKERS parallel R processes on one node.
# Each worker processes every Nth gene (stride pattern).
# After all workers finish, merges results into the final output.
#
# Usage: GWAS_NAME=UKBB_ALT N_WORKERS=4 sbatch --array=1 src/06_susie_coloc_bigmem.sh
#        (or --array=1-2,5,8,18 for multiple chromosomes)

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

export CHR_FILTER="${SLURM_ARRAY_TASK_ID}"
N_WORKERS="${N_WORKERS:-4}"

echo "=== SuSiE-COLOC Bigmem: ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "Workers: ${N_WORKERS} | CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Start: $(date)"

# Launch N workers in parallel
PIDS=()
for w in $(seq 0 $((N_WORKERS - 1))); do
  echo "  Launching worker ${w}/${N_WORKERS}..."
  WORKER_ID=${w} N_WORKERS=${N_WORKERS} \
    Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}" &
  PIDS+=($!)
done

echo "All ${N_WORKERS} workers launched. PIDs: ${PIDS[*]}"

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
WORKER_FILES=""
for w in $(seq 0 $((N_WORKERS - 1))); do
  wf="${OUT_DIR}/worker_chr${CHR_FILTER}_w${w}.csv"
  if [ -f "$wf" ]; then
    WORKER_FILES="${WORKER_FILES} ${wf}"
  fi
done

if [ -n "$WORKER_FILES" ]; then
  echo "Merging worker outputs..."
  # Use R to merge and deduplicate
  Rscript -e "
    library(data.table)
    files <- commandArgs(trailingOnly = TRUE)
    dt_list <- lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL))
    dt_list <- dt_list[!sapply(dt_list, is.null)]
    if (length(dt_list) > 0) {
      merged <- rbindlist(dt_list, fill = TRUE)
      merged <- merged[!duplicated(ensembl)]
      fwrite(merged, '${OUT_FILE}')
      cat('Merged', nrow(merged), 'genes into ${OUT_FILE}\n')
    }
  " ${WORKER_FILES}

  # Clean up worker files
  rm -f ${WORKER_FILES}
  rm -f ${OUT_DIR}/worker_chr${CHR_FILTER}_w*.csv
  rm -f ${OUT_DIR}/checkpoint_chr${CHR_FILTER}_w*.csv
fi

# Clean up checkpoint
rm -f "${OUT_DIR}/checkpoint_chr${CHR_FILTER}.csv"

echo "End: $(date)"
