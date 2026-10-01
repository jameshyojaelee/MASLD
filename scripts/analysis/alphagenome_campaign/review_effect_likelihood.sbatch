#!/bin/bash
#SBATCH --job-name=masld_review_likelihood
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --dependency=afterok:21774086
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_likelihood_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_likelihood_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_effect_likelihood.py" \
 --run "$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/effect_likelihood_21774086" \
 --out "$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/effect-likelihood-$SLURM_JOB_ID" --draws 2000
