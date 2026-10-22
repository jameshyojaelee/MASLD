#!/bin/bash
#SBATCH --job-name=H36_stage60
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-codepath_test_stage60-%j.out

# Code-path test of the whole of h36_03_stage_test.R (including the calibration and cohort-check tail) on the first
# 60 clusters. Not a result.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
H=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
test ! -e "$H/codepath_test_stage60"; mkdir -p "$H/codepath_test_stage60"
$RS $SC/h36_03_stage_test.R $REPO $H/build $H/codepath_test_stage60 10 60
