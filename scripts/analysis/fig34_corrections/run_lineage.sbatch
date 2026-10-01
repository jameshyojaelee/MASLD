#!/usr/bin/env bash
#SBATCH --job-name=fig34-lineage
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python scripts/analysis/fig34_corrections/03_lineage.py "$1/lineage_matched"
