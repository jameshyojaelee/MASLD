#!/bin/bash
#SBATCH --job-name=u15_s2f
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/logs/%j.out
set -euo pipefail
U=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $U/out_s2ai_carma
$PY -m pip freeze > $U/out_s2ai_carma/pip_freeze.txt
sha256sum $U/u15_s2a_i.py $U/u15_s2c.py > $U/out_s2ai_carma/code_sha256.txt
$PY -B $U/u15_s2a_i.py $U/out_s2ai_carma carma
echo DONE $(date -u)
