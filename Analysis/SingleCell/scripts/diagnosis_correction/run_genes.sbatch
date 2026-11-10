#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/genes_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/genes_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f "$TASK/PREPARED"
"$RS" Analysis/SingleCell/scripts/pseudobulk_de.R
touch "$TASK/BINARY_DE_COMPLETED"
export SINGLECELL_BINARY_DE_OUT="$TASK/pseudobulk_de_strict_control"
export SINGLECELL_STRICT_CONTROL=1
"$RS" Analysis/SingleCell/scripts/pseudobulk_de.R
touch "$TASK/STRICT_BINARY_DE_COMPLETED"
