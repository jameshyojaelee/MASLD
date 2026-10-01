#!/bin/bash
#SBATCH --job-name=masld_f5g_gc_sensitivity
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5g_gc_sensitivity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5g_gc_sensitivity_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
DATA=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data
cd "$PROJ/scripts/analysis/alphagenome_campaign"
"$PY" data_p4_f5_gc_sensitivity.py \
  --reduced "$DATA/p4-full-21772842" "$DATA/u2c-reduce-21783091" "$DATA/f5background-reduce-21783090" \
  --out "$DATA/f5g-gc-sensitivity-${SLURM_JOB_ID}"
