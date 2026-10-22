#!/bin/bash
#SBATCH --job-name=matrix_rg_pairs
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=logs/matrix_rg_pairs_%j.out
#SBATCH --error=logs/matrix_rg_pairs_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
MATRIX_RUN_ROOT="GWAS/ldsc/rg_missing_compatible_$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
mkdir -p "$MATRIX_RUN_ROOT"
echo "$MATRIX_RUN_ROOT"
cp scripts/analysis/run_ldsc_matrix_missing_pairs.sbatch "$MATRIX_RUN_ROOT/"
cp scripts/analysis/complete_ldsc_matrix_pairs.py "$MATRIX_RUN_ROOT/"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/complete_ldsc_matrix_pairs.py "$MATRIX_RUN_ROOT"
