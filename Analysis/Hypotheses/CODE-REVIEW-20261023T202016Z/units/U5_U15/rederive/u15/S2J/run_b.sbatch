#!/bin/bash
#SBATCH --job-name=u15_b
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2J/s2j_s2m_rederive.py 2>&1 | tee /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2J/s2j_s2m_rederive.log
