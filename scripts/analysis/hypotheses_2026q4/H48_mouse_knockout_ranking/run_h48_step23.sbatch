#!/bin/bash
#SBATCH --job-name=H48_fit_score
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=%x_%j.log
set -euo pipefail
OUT=$1
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H48_mouse_knockout_ranking
PY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_gears/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
$RS $S/h48_step2_fit.R $OUT 16
$PY -I $S/h48_step3_score.py --out $OUT
