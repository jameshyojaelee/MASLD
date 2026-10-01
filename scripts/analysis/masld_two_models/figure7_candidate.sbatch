#!/bin/bash
#SBATCH --job-name=masld_two_models_figure7
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:21849197:21849200
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_two_models_figure7_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_two_models_figure7_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$PROJ/scripts/analysis/masld_two_models/figure7_candidate.py" \
  --chrom "$BASE/chromatin_conditional_21849192" \
  --acquisition "$BASE/chromatin_acquisition_21849197" \
  --allele "$BASE/allele_combination_21849200" \
  --out "$PROJ/figures/candidates/figure7_two_models_20260922_${SLURM_JOB_ID}"
