#!/bin/bash
#SBATCH --job-name=masld_allele_combination
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:21849187
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_allele_combination_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_allele_combination_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$PROJ/scripts/analysis/masld_two_models/allele_combination.py" \
  --inputs "$BASE/allele_stack_inputs_21849187" \
  --out "$BASE/allele_combination_${SLURM_JOB_ID}"
