#!/bin/bash
#SBATCH --job-name=x1-selection-fix
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/xenium/logs/%x_%j.out
#SBATCH --error=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/xenium/logs/%x_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
P=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
OUT=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/xenium/results_${SLURM_JOB_ID}
mkdir "$OUT"
"$P" -I Analysis/Hypotheses/X1-20261028T103410Z/scripts/x1_xenium_v2.py "$OUT"
"$P" -m pip freeze > "$OUT/pip_freeze.txt"
