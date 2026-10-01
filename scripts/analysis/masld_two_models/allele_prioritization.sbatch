#!/bin/bash
#SBATCH --job-name=masld_allele_prioritization
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:21849200
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_allele_prioritization_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_allele_prioritization_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$PROJ/scripts/analysis/masld_two_models/allele_prioritization.py" \
  --predictions "$BASE/allele_combination_21849200/matched_development_predictions.tsv.gz" \
  --out "$BASE/allele_prioritization_${SLURM_JOB_ID}"
