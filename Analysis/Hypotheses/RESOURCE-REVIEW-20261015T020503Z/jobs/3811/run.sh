#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3811
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
mkdir -p $TMPDIR
cd $J
date
JOBDIR=$J /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $J/feasibility.R
date
