#!/bin/bash
#SBATCH --job-name=qualified_CCC
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=logs/qualified_CCC_%j.out
#SBATCH --error=logs/qualified_CCC_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export MASLD_PROJECT_ROOT="$PWD"
export CCC_QUALIFIED_ROOT="$PWD/Analysis/SingleCell/candidates/source-qualified-CCC-20261014_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/analysis/bulk_refit_followup/refit_source_qualified_ccc.R
cp scripts/analysis/bulk_refit_followup/refit_source_qualified_ccc.R scripts/analysis/bulk_refit_followup/run_source_qualified_ccc.sbatch "$CCC_QUALIFIED_ROOT/"
