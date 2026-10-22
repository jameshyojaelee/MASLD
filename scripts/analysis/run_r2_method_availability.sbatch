#!/bin/bash
#SBATCH --job-name=r2_availability
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=logs/r2_availability_%j.out
#SBATCH --error=logs/r2_availability_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export MASLD_PROJECT_ROOT="$PWD"
export R2_AVAILABILITY_ROOT="$PWD/GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/method_availability_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/audit_r2_method_availability.py
cp scripts/analysis/audit_r2_method_availability.py scripts/analysis/run_r2_method_availability.sbatch "$R2_AVAILABILITY_ROOT/"
