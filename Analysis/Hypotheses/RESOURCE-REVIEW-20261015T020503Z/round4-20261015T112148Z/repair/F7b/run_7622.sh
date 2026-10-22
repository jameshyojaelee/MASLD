#!/bin/bash
#SBATCH --job-name=7622
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b/logs/%x-%j.err
# Re-render 7D from the 7621 table; only change: caption prints tie-ceiling ratios at 3 decimals
# (CVH 0.680, not 0.6805, so it cannot be read as 0.681).
set -euo pipefail
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b
FIGOUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair
FIG7_SCRIPTS="$W" FIG7_EXEC_DIR="$FIGOUT" FIG7_OUT_DIR="$FIGOUT" \
  "$ENV/Rscript" "$W/fig7d_measured_bounds_and_join_structure.R" 2>&1 | tee "$W/logs/fig7d_render_7622.out"
grep '^\[CAPTION\]' "$W/logs/fig7d_render_7622.out" > "$FIGOUT/fig7D_caption.txt"
pdftotext -layout "$FIGOUT/fig7D_measured_bounds_and_join_structure.pdf" "$W/logs/fig7d_pdftotext_7622.txt" || true
pdftoppm -r 220 -png "$FIGOUT/fig7D_measured_bounds_and_join_structure.pdf" "$W/logs/fig7d_check_7622" || true
