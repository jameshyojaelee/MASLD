#!/bin/bash
#SBATCH --job-name=H35_test_bb
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H35-20261009T180343Z/logs/test_bb-%j.out
#SBATCH --error=Analysis/Hypotheses/H35-20261009T180343Z/logs/test_bb-%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/analysis/hypotheses_2026q4/H35_editing/test_lib_bb.R
