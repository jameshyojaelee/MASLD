#!/bin/bash
#SBATCH --job-name=masld_review_RNA_masks
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_RNA_masks_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_RNA_masks_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
 "$review_root/scripts/analysis/alphagenome_campaign/review_donor_mask_precheck.py" \
 --out "$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/donor-mask-precheck-$SLURM_JOB_ID"
