#!/bin/bash
#SBATCH --job-name=masld_review_f6
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_f6_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_f6_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
WK=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
cd "$PROJ/scripts/analysis/alphagenome_campaign"
"$PY" review_f6_overlap.py \
  --reduced "$WK/data/p4-full-21772842" "$WK/data/u2c-reduce-21783091" "$WK/data/f5background-reduce-21783090" \
  --producer "$WK/data/f5-f6-families-21797329" \
  --out "$WK/review/f6-overlap-${SLURM_JOB_ID}"
