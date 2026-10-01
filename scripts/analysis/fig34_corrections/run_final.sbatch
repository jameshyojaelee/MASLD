#!/usr/bin/env bash
#SBATCH --job-name=fig34-final
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$PY" scripts/analysis/fig34_corrections/05_figures_transport.py "$1"
"$PY" scripts/analysis/fig34_corrections/08_verify_export.py "$1"
