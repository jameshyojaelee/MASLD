#!/bin/bash
#SBATCH --job-name=u8a_qc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
U=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I $U/07_qc_donor_scores.py $U/out/qc_donor_scores
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > $U/out/qc_donor_scores/pip_freeze.txt
