#!/bin/bash
#SBATCH --job-name=ag2x_c00_cellsets
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.err
set -euo pipefail
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export PYTHONHASHSEED=20261111 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
mkdir -p $O/stageC/inputs/cellsets
$PY -m pip freeze > $O/stageC/inputs/pip_freeze_c00_rapids_singlecell.txt
$PY -I $O/code/c00_build_cellsets.py
sha256sum $O/code/c00_build_cellsets.py $O/stageC/inputs/cellsets/*.h5ad $O/stageC/inputs/variants_stageC.tsv > $O/stageC/inputs/c00_sha256.txt
