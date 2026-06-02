#!/bin/bash -l
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=72:00:00
#SBATCH --array=1-22
#SBATCH --job-name=coloc_pukbb_sex
#SBATCH --output=GWAS/finemapping/logs/coloc_pukbb_sex_%x_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/coloc_pukbb_sex_%x_%A_%a.err

# ---------------------------------------------------------------------------
# Per-chromosome SuSiE-COLOC for Pan-UKBB/Neale sex-stratified GWAS
# (Team B1, headline analysis A1).
#
# Set externally:  GWAS_NAME ∈ {PanUKBB_F_ALT, PanUKBB_F_AST, PanUKBB_F_GGT,
#                                PanUKBB_M_ALT, PanUKBB_M_AST, PanUKBB_M_GGT}
# LD reference:    PolyFun EUR (default LD_PANEL=polyfun, canonical 2026-05-06)
# eQTL ref:        Broadaway liver (canonical, EUR N=1,183)
# ---------------------------------------------------------------------------

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
# 06_susie_coloc.sh upstream uses `finemapping`; honour the canonical pipeline env.
micromamba activate finemapping

export CHR_FILTER="${SLURM_ARRAY_TASK_ID}"
export LD_PANEL="${LD_PANEL:-polyfun}"
export N_WORKERS="${N_WORKERS:-1}"

if [ -z "${GWAS_NAME:-}" ]; then
  echo "ERROR: GWAS_NAME env var must be set (e.g. PanUKBB_F_ALT)"
  exit 2
fi

echo "=== SuSiE-COLOC sex-stratified: ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "  LD panel: ${LD_PANEL}"
echo "  N_WORKERS: ${N_WORKERS}"
echo "  Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Start: $(date)"

if [ "${N_WORKERS}" -le 1 ]; then
  Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}"
else
  PIDS=()
  for w in $(seq 0 $((N_WORKERS - 1))); do
    WORKER_ID=${w} N_WORKERS=${N_WORKERS} \
      Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}" &
    PIDS+=($!)
  done
  echo "  Launched ${N_WORKERS} workers. PIDs: ${PIDS[*]}"

  FAIL=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=$((FAIL + 1))
  done
  if [ $FAIL -eq $N_WORKERS ]; then
    echo "ERROR: all workers failed"; exit 1
  fi

  OUT_DIR="results/susie_coloc/${GWAS_NAME}"
  OUT_FILE="${OUT_DIR}/susie_coloc_chr${CHR_FILTER}.csv"
  WORKER_FILES=$(ls ${OUT_DIR}/worker_chr${CHR_FILTER}_w*.csv 2>/dev/null || true)
  if [ -n "$WORKER_FILES" ]; then
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
    rm -f ${OUT_DIR}/worker_chr${CHR_FILTER}_w*.csv 2>/dev/null || true
    rm -f ${OUT_DIR}/checkpoint_chr${CHR_FILTER}_w*.csv 2>/dev/null || true
  fi
fi

rm -f "results/susie_coloc/${GWAS_NAME}/checkpoint_chr${CHR_FILTER}.csv" 2>/dev/null || true

echo "End: $(date)"
