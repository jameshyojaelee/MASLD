#!/usr/bin/env bash
#SBATCH --job-name=fig34-tables
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
FIG34_ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
"$FIG34_ENV/bin/python" scripts/analysis/fig34_corrections/01_continuum.py "$1/continuum"
"$FIG34_ENV/bin/Rscript" scripts/analysis/fig34_corrections/02_protein.R "$1/protein"
