#!/bin/bash
#SBATCH --job-name=masld_review_zero_outputs
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --dependency=afterok:21773146
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_zero_outputs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_zero_outputs_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
zero_review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
zero_review_base=$zero_review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
zero_review_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$zero_review_python" "$zero_review_root/scripts/analysis/alphagenome_campaign/review_zero_outputs.py" \
 --audit "$zero_review_base/model/zero_outputs_21773146" \
 --out "$zero_review_base/review/zero-outputs-$SLURM_JOB_ID"
