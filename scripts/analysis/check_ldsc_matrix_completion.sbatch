#!/bin/bash
#SBATCH --job-name=rg_matrix_check
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=logs/rg_matrix_check_%j.out
#SBATCH --error=logs/rg_matrix_check_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
MATRIX_CHECK_ROOT="$1"
cp scripts/analysis/check_ldsc_matrix_completion.py "$MATRIX_CHECK_ROOT/check_ldsc_matrix_completion.py"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/check_ldsc_matrix_completion.py "$MATRIX_CHECK_ROOT"
