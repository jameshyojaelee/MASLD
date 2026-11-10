#!/bin/bash
#SBATCH --job-name=eqtl_dedup_LD
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-20261110T023251Z/logs/ld_diag_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-20261110T023251Z/logs/ld_diag_%j.err
set -eo pipefail
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
export RUN_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-20261110T023251Z
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design LD_PANEL=polyfun
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
Rscript ${RUN_DIR}/code/ld_diagnostic.R
