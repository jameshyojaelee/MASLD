#!/bin/bash
#SBATCH --job-name=H34_y_loss_r2
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H34-20261008T015507Z/slurm-%j.out

set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H34-20261008T015507Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=34 OMP_NUM_THREADS=4
test -f "$OUT/PRESPEC_H34.md"
sha256sum "$OUT/PRESPEC_H34.md" "$OUT/DEVIATIONS_H34.md" "$REPO/scripts/analysis/hypotheses_2026q4/H34_y_loss/h34_y_loss.py"
"$PY" -m pip freeze > "$OUT/pip_freeze.txt"
"$PY" "$REPO/scripts/analysis/hypotheses_2026q4/H34_y_loss/h34_y_loss.py" "$REPO" "$OUT" "$REPO/Analysis/Hypotheses/H34-20261008T014810Z/intermediate"
