#!/bin/bash -l
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=1900G
#SBATCH --time=90:00:00  # capped at 90h per CLAUDE.md absolute-max policy (T2.3 2026-04-22); per-chr resume handles re-runs
#SBATCH --job-name=susie_coloc
#SBATCH --output=logs/susie_coloc_bigmem_multi_%j.out
#SBATCH --error=logs/susie_coloc_bigmem_multi_%j.err

# Bigmem multi-chr SuSiE-COLOC: processes MULTIPLE chromosomes sequentially,
# each with N_WORKERS parallel R processes. Maximizes bigmem node utilization.
#
# Usage: GWAS_NAME=UKBB_ALT CHR_LIST="1,2,5,8,18" N_WORKERS=16 sbatch src/06_susie_coloc_bigmem_multi.sh

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

N_WORKERS="${N_WORKERS:-16}"
CHR_LIST="${CHR_LIST:-1,2,5,8,18}"

# Validate GWAS_NAME
if [ -z "${GWAS_NAME}" ]; then
  echo "ERROR: GWAS_NAME not set"
  exit 1
fi

echo "=== SuSiE-COLOC Bigmem Multi-Chr: ${GWAS_NAME} ==="
echo "Chromosomes: ${CHR_LIST}"
echo "Workers per chr: ${N_WORKERS}"
echo "Job: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo ""

# Process each chromosome sequentially
IFS=',' read -ra CHRS <<< "$CHR_LIST"
for chr in "${CHRS[@]}"; do
  export CHR_FILTER="${chr}"
  OUT_DIR="results/susie_coloc/${GWAS_NAME}"
  OUT_FILE="${OUT_DIR}/susie_coloc_chr${chr}.csv"

  # Skip if already completed (has SuSiE columns)
  if [ -f "$OUT_FILE" ]; then
    has_susie=$(head -1 "$OUT_FILE" | grep -c "PP.H4.susie" || true)
    if [ $has_susie -gt 0 ]; then
      echo "=== chr${chr}: SKIPPING (already has SuSiE results) ==="
      continue
    fi
  fi

  echo "=== chr${chr}: Starting (${N_WORKERS} workers) ==="
  echo "  $(date)"

  # Launch N workers in parallel
  PIDS=()
  for w in $(seq 0 $((N_WORKERS - 1))); do
    WORKER_ID=${w} N_WORKERS=${N_WORKERS} \
      Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${chr}" &
    PIDS+=($!)
  done

  # Wait for all workers
  FAIL=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=$((FAIL + 1))
  done

  if [ $FAIL -eq $N_WORKERS ]; then
    echo "  ERROR: ALL ${N_WORKERS} workers failed for chr${chr}"
    continue
  elif [ $FAIL -gt 0 ]; then
    echo "  WARNING: ${FAIL}/${N_WORKERS} workers failed for chr${chr}"
  fi

  # Clean stale worker files from any prior failed run BEFORE merging
  # (prevents contamination from old partial results)

  # Merge worker outputs
  WORKER_FILES=$(ls ${OUT_DIR}/worker_chr${chr}_w*.csv 2>/dev/null || true)
  if [ -n "$WORKER_FILES" ]; then
    echo "  Merging $(echo "$WORKER_FILES" | wc -w) worker files..."
    Rscript -e "
      library(data.table)
      files <- commandArgs(trailingOnly = TRUE)
      dt_list <- lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL))
      dt_list <- dt_list[!sapply(dt_list, is.null)]
      if (length(dt_list) > 0) {
        merged <- rbindlist(dt_list, fill = TRUE)
        merged <- merged[!duplicated(ensembl)]
        out_file <- '${OUT_FILE}'
        fwrite(merged, out_file)
        cat('  Merged', nrow(merged), 'genes into', out_file, '\n')
      } else {
        cat('  WARNING: No valid worker files to merge\n')
      }
    " ${WORKER_FILES}
    rm -f ${OUT_DIR}/worker_chr${chr}_w*.csv 2>/dev/null || true
    rm -f ${OUT_DIR}/checkpoint_chr${chr}_w*.csv 2>/dev/null || true
  else
    echo "  WARNING: No worker output files found for chr${chr}"
  fi
  rm -f "${OUT_DIR}/checkpoint_chr${chr}.csv" 2>/dev/null || true

  echo "  chr${chr}: DONE at $(date)"
  echo ""
done

echo "=== All chromosomes complete ==="
echo "End: $(date)"
