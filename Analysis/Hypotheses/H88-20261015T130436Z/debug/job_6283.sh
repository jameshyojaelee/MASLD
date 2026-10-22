#!/bin/bash
#SBATCH --job-name=6283
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H88-20261015T130436Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H88-20261015T130436Z/logs/%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H88-20261015T130436Z
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$RUN/tmp XDG_CACHE_HOME=$RUN/tmp MPLCONFIGDIR=$RUN/tmp OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
mkdir -p $TMPDIR
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $RUN/run_h88.py $RUN/debug/d3_assert_check --debug --bscale 0.02 --cache-in $RUN/tmp/pass1_cache.pkl
