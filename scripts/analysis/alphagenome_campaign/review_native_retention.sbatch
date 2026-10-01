#!/bin/bash
#SBATCH --job-name=masld_review_native_ATAC
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:21774009
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_native_ATAC_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_native_ATAC_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
review_base=$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_native_retention.py" \
 --fixture "$review_base/model/native_retention_panel_21773996" \
 --run "$review_base/model/native_retention_21774009" \
 --out "$review_base/review/native-retention-$SLURM_JOB_ID"
