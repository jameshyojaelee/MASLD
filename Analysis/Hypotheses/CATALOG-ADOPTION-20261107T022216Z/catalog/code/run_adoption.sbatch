#!/bin/bash
#SBATCH --job-name=masld-catalog-adopt
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/catalog/logs/%j.out
#SBATCH --error=Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/catalog/logs/%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
TASK=$PWD/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/catalog
"$PY" scripts/portal/build_adopted_eligibility_catalog.py --project-root "$PWD"
"$PY" -m pip freeze > "$TASK/logs/python_environment.txt"
