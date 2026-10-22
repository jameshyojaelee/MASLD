#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4071
mkdir -p $J/tmp
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp
cd $J
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript check_hk.R
