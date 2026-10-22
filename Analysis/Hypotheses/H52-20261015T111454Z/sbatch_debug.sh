#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H52-20261015T111454Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H52-20261015T111454Z/logs/%j.err
# Debug run: H52_DEBUG=1 permutes every GSE213621 value per gene before use; outputs to tmp/debug_<tag>.
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H52-20261015T111454Z
export TMPDIR=$RUN/tmp XDG_CACHE_HOME=$RUN/tmp MPLCONFIGDIR=$RUN/tmp OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
export H52_DEBUG=1
mkdir -p $TMPDIR
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $RUN/run_h52.R
