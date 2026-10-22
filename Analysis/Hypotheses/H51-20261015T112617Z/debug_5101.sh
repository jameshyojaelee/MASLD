#!/bin/bash
#SBATCH --job-name=5101
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z/logs/%j.err
set -euo pipefail
H=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR=$H/tmp XDG_CACHE_HOME=$H/tmp MPLCONFIGDIR=$H/tmp
export H51_DEBUG=1
cd $H
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript stage1_blinded.R
S1=$(ls -d $H/debug-*/stage1 | tail -1)
echo "stage1 dir $S1"
H51_STAGE1_DIR=$S1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript stage2_outcome.R
