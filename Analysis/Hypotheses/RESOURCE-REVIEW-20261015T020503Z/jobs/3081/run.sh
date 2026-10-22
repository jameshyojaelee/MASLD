#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3081
export TMPDIR=$J/tmp; mkdir -p $TMPDIR
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
cd $J
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript pc_check.R
