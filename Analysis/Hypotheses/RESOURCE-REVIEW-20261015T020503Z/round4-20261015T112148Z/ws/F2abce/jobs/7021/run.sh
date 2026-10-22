#!/bin/bash
set -euo pipefail
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F2abce
export TMPDIR=$W/tmp XDG_CACHE_HOME=$W/tmp
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $W/scripts/fig2b_same_row.R > $W/jobs/7021/fig2b.log 2>&1
