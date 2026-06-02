#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --job-name=MASLD-COLOC
#SBATCH --output=logs/orphan_array_%A_%a.out
#SBATCH --error=logs/orphan_array_%A_%a.err
#SBATCH --array=1-25%12

# CPU array variant of 06_susie_coloc_bigmem_recovery.sh
# Each array task processes ONE (GWAS, CHR) pair with N_WORKERS=8
# Reads line $SLURM_ARRAY_TASK_ID from TASK_FILE (format: GWAS|CHR)
#
# Usage:
#   TASK_FILE=data/recovery/tasks_orphan_array.txt sbatch src/06_susie_coloc_cpu_array.sh

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/finemapping/bin/Rscript"
N_WORKERS="${N_WORKERS:-8}"
TASK_FILE="${TASK_FILE:-data/recovery/tasks_orphan_array.txt}"

if [ ! -f "$TASK_FILE" ]; then
  echo "ERROR: TASK_FILE not found: $TASK_FILE"
  exit 1
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_FILE")
GWAS=$(echo "$LINE" | cut -d'|' -f1)
CHR=$(echo "$LINE" | cut -d'|' -f2)

if [ -z "$GWAS" ] || [ -z "$CHR" ]; then
  echo "ERROR: Could not parse line ${SLURM_ARRAY_TASK_ID}: '$LINE'"
  exit 1
fi

OUT_DIR="results/susie_coloc/${GWAS}"
OUT_FILE="${OUT_DIR}/susie_coloc_chr${CHR}.csv"
mkdir -p "$OUT_DIR"

echo "=== Array task ${SLURM_ARRAY_TASK_ID}: ${GWAS} chr${CHR} ==="
echo "Workers: ${N_WORKERS}"
echo "Mem: ${SLURM_MEM_PER_NODE}M"
echo "Start: $(date)"

# Skip if already SuSiE-complete
if [ -f "$OUT_FILE" ] && head -1 "$OUT_FILE" | grep -q "PP.H4.susie"; then
  echo "Already SuSiE-complete — exiting"
  exit 0
fi

# Clean stale worker files from prior runs with different N_WORKERS
existing=$(ls "${OUT_DIR}/worker_chr${CHR}_w"*.csv 2>/dev/null | wc -l || echo 0)
if [ "$existing" -gt 0 ] && [ "$existing" -ne "$N_WORKERS" ]; then
  echo "Removing ${existing} stale worker files (mismatch vs N_WORKERS=${N_WORKERS})"
  rm -f "${OUT_DIR}/worker_chr${CHR}_w"*.csv
  rm -f "${OUT_DIR}/checkpoint_chr${CHR}_w"*.csv
  rm -f "${OUT_DIR}/checkpoint_chr${CHR}.csv"
fi

# Launch N_WORKERS parallel R processes
PIDS=()
for w in $(seq 0 $((N_WORKERS - 1))); do
  WORKER_ID=${w} N_WORKERS=${N_WORKERS} CHR_FILTER=${CHR} \
    ${RSCRIPT} src/06_susie_coloc.R "${GWAS}" "${CHR}" \
    > "logs/orphan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_w${w}.out" 2>&1 &
  PIDS+=($!)
done

FAIL=0
for pid in "${PIDS[@]}"; do
  wait "$pid" || FAIL=$((FAIL + 1))
done

echo "Workers finished. Failures: ${FAIL}/${N_WORKERS}"

# Merge worker outputs
WORKER_FILES=$(ls ${OUT_DIR}/worker_chr${CHR}_w*.csv 2>/dev/null || true)
if [ -n "$WORKER_FILES" ]; then
  echo "Merging $(echo "$WORKER_FILES" | wc -w) worker files..."
  ${RSCRIPT} -e "
    suppressPackageStartupMessages(library(data.table))
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
fi

# Verify & cleanup
if head -1 "$OUT_FILE" 2>/dev/null | grep -q "PP.H4.susie"; then
  rm -f "${OUT_DIR}/worker_chr${CHR}_w"*.csv
  rm -f "${OUT_DIR}/checkpoint_chr${CHR}_w"*.csv
  rm -f "${OUT_DIR}/checkpoint_chr${CHR}.csv"
  echo "${GWAS} chr${CHR}: DONE at $(date)"
else
  echo "${GWAS} chr${CHR}: merge missing PP.H4.susie — worker files preserved for retry"
  exit 1
fi
