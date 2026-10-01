#!/bin/bash
#SBATCH --job-name=masld_p4_background
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_p4_background_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_p4_background_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
background_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$background_root/scripts/analysis/alphagenome_campaign/data_p4_background.py" \
  --out "$background_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/p4-background-$SLURM_JOB_ID"
