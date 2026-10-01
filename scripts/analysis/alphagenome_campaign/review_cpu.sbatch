#!/bin/bash
#SBATCH --job-name=ag_week1_b2_review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/job-%j.out
#SBATCH --error=GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/job-%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONNOUSERSITE=1
REVIEW_PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
REVIEW_OUT=GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review
"$REVIEW_PY" scripts/analysis/alphagenome_campaign/review_evidence.py --output "$REVIEW_OUT/evidence-${SLURM_JOB_ID}"
"$REVIEW_PY" scripts/analysis/alphagenome_campaign/review_donor.py --output "$REVIEW_OUT/b2-${SLURM_JOB_ID}"
