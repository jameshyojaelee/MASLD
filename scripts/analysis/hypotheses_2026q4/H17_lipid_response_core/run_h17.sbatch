#!/bin/bash
#SBATCH --job-name=H17_core
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
# H17: lipid-response core shared by iPSC-hepatocytes (GSE342568) and in vivo
# hepatocyte donor pseudobulk. Single job, no array.
# Usage: sbatch --output=<HYP_DIR>/logs/%x_%j.out --error=<HYP_DIR>/logs/%x_%j.err run_h17.sbatch <HYP_DIR>
set -euo pipefail
HYP="${1:?hypothesis dir required}"
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H17_lipid_response_core
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PATH="$ENV/bin:$PATH"
test -f "$HYP/PRESPEC_H17.md"
OUT="$HYP/results_$(date -u +%Y%m%dT%H%M%SZ)_job${SLURM_JOB_ID}"
test ! -e "$OUT"
mkdir -p "$OUT"
echo "OUT=$OUT"
cd "$OUT"
"$ENV/bin/Rscript" "$HERE/01_de_and_concordance.R" "$HYP" "$OUT"
if ls "$OUT"/STOP_* >/dev/null 2>&1; then
  echo "stopped by prespecified rule"; cat "$OUT"/STOP_*
  "$ENV/bin/Rscript" -e 'suppressPackageStartupMessages({library(limma);library(edgeR);library(data.table)}); capture.output(sessionInfo(), file=commandArgs(TRUE)[1])' "$OUT/sessionInfo.txt"
  exit 0
fi
"$ENV/bin/Rscript" "$HERE/02_resampling.R" "$HYP" "$OUT"
"$ENV/bin/Rscript" "$HERE/03_sensitivity_and_call.R" "$HYP" "$OUT"
sha256sum "$HERE"/*.R "$HERE"/run_h17.sbatch > "$OUT/script_sha256.txt"
echo "H17 done"
