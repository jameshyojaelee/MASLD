#!/bin/bash
#SBATCH --job-name=WP5_obs
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z/logs/%j.err
# WP5: 117-program observability panel (detection rates only; candidate).
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=8
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
sha256sum "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z/scripts/01_observability_panel.py" > "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z/script_sha256.txt"
"$PY" -I "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z/scripts/01_observability_panel.py"
"$PY" -I -m pip freeze > "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z/pip_freeze_spatial.txt"
