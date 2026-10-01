#!/bin/bash
#SBATCH --job-name=masld_review_adapt_precision
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_adapt_precision_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_adapt_precision_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_base=$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
review_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_adaptation_precision.py" \
 --sweep "$review_base/model/sweep_21772726" \
 --manifest "$review_base/model/native/sequence_manifest.tsv" \
 --out "$review_base/review/adaptation-precision-$SLURM_JOB_ID" --draws 2000
