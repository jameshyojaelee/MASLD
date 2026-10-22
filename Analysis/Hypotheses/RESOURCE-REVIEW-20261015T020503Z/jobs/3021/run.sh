#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3021
mkdir -p $J/tmp
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp OMP_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $J/null_power.R > $J/null_power.log 2>&1
echo done
