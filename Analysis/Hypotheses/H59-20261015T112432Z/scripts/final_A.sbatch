#!/bin/bash
#SBATCH --job-name=5900
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H59-20261015T112432Z/logs/%j.out
#SBATCH --error=Analysis/Hypotheses/H59-20261015T112432Z/logs/%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H59-20261015T112432Z
export TMPDIR=$R/tmp H59_OUT=$R/out
bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
"$bin/Rscript" --vanilla $R/scripts/01_h59_gates_null.R
