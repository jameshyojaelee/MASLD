#!/bin/bash
#SBATCH --job-name=cam_step
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=1
cd "$MASLD_PROJECT_ROOT"
RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PYSPATIAL=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
case "$1" in
  *.R)  "$RSCRIPT" --vanilla "$1" ;;
  *.py) "$PYSPATIAL" "$1" --out-root "$CAM_OUT_ROOT" ;;
esac
