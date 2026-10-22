#!/bin/bash
#SBATCH --job-name=verify_gse126848
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/verify_gse126848_%j.out
#SBATCH --error=logs/verify_gse126848_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export BULK_FOLLOWUP_ROOT="$PWD/figures/candidates/bulk-context-774-20261014T140904Z_22236311"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cp scripts/analysis/bulk_refit_followup/verify_native_gse126848.R "$BULK_FOLLOWUP_ROOT/source/"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla "$BULK_FOLLOWUP_ROOT/source/verify_native_gse126848.R"
