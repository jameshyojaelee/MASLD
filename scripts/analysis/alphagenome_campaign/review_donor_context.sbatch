#!/bin/bash
#SBATCH --job-name=masld_review_donor_context
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:21773708
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_donor_context_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_donor_context_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_base=$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
 "$review_root/scripts/analysis/alphagenome_campaign/review_donor_context.py" \
 --run "$review_base/model/donor_context_21773708" \
 --out "$review_base/review/donor-context-$SLURM_JOB_ID"
