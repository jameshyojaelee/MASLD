#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4201
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp
cd $J
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript verify_g2_02.R
