#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/coloc_coding_noncoding/logs/render_03_05_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/coloc_coding_noncoding/logs/render_03_05_%j.err
# Re-render coloc coding/non-coding chain 03->04->05 on C2 canonical.
# 01/02 intermediates already exist (regenerated 2026-06-08 12:42-43); only
# 03 (aggregate + DEG overlap) -> 04 (F1-F8 panels) -> 05 (fig3 main) needed.
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -uo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=$BASE/scripts/figures/coloc_coding_noncoding
export MASLD_PROJECT_ROOT="$BASE"
mkdir -p $SCRIPTS/logs
cd $BASE

rc=0
for s in 03_aggregate_and_deg_overlap.R 04_render_figures.R 05_render_fig3_main.R; do
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running $s ..."
  echo "=========================================="
  if Rscript "$SCRIPTS/$s"; then
    echo "[$(date '+%H:%M:%S')] $s SUCCEEDED"
  else
    echo "[$(date '+%H:%M:%S')] $s FAILED"
    rc=1
  fi
done
echo "[$(date)] DONE rc=$rc"
exit $rc
