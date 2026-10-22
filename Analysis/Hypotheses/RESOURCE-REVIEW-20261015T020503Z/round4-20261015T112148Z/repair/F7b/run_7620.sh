#!/bin/bash
#SBATCH --job-name=7620
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b/logs/%x-%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b
FIGOUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair
cd "$W"
echo "host $(hostname) start $(date -Is)"
"$ENV/python" r0802_transfer_primary_units.py 2>&1 | tee "$W/logs/r0802_transfer_primary_units.out"
FIG7_SCRIPTS="$W" FIG7_EXEC_DIR="$FIGOUT" FIG7_OUT_DIR="$FIGOUT" \
  "$ENV/Rscript" "$W/fig7d_measured_bounds_and_join_structure.R" 2>&1 | tee "$W/logs/fig7d_render.out"
grep '^\[CAPTION\]' "$W/logs/fig7d_render.out" > "$FIGOUT/fig7D_caption.txt" || true
pdftotext -layout "$FIGOUT/fig7D_measured_bounds_and_join_structure.pdf" "$W/logs/fig7d_pdftotext.txt" || true
pdftoppm -r 200 -png "$FIGOUT/fig7D_measured_bounds_and_join_structure.pdf" "$W/logs/fig7d_check" || true
"$ENV/python" r1101_separability_forms_sim.py "$FIGOUT" 2>&1 | tee "$W/logs/r1101_forms_sim.out"
"$ENV/python" -m pip freeze > "$FIGOUT/pip_freeze.txt" 2>/dev/null || true
SI="$FIGOUT/sessionInfo.txt" "$ENV/Rscript" -e 'suppressPackageStartupMessages({library(ggplot2);library(patchwork)}); writeLines(capture.output(sessionInfo()), Sys.getenv("SI"))'
echo "end $(date -Is)"
