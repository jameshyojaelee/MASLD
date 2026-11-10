#!/bin/bash
#SBATCH --job-name=unique-refit
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/NUCLEUS-FIX-20261107/refit_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MASLD_PROJECT_ROOT="$PWD"
task="$PWD/Analysis/Hypotheses/NUCLEUS-FIX-20261107"
export SINGLECELL_ATLAS_H5AD="$task/unique_counts/unique_counts.h5ad"
export SINGLECELL_DIAGNOSIS_ROOT="$task/unique_counts/metadata"
export SINGLECELL_PSEUDOBULK_DIR="$task/unique_counts/pseudobulk"
export SINGLECELL_BINARY_DE_OUT="$task/binary"
export SINGLECELL_STAGE_DE_OUT="$task/stage"
export SINGLECELL_ALLCELL_OUT="$task/allcell"
rs=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
"$rs" Analysis/SingleCell/scripts/pseudobulk_de.R
"$rs" Analysis/SingleCell/scripts/run_pseudobulk_de_coarse_stages.R
"$rs" Analysis/SingleCell/scripts/314_allcell_pseudobulk_de.R
export SINGLECELL_BINARY_DE_OUT="$task/binary_strict"
export SINGLECELL_STRICT_CONTROL=1
"$rs" Analysis/SingleCell/scripts/pseudobulk_de.R
