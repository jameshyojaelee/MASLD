#!/bin/bash
set -euo pipefail
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F2abce
export TMPDIR=$W/tmp XDG_CACHE_HOME=$W/tmp MPLCONFIGDIR=$W/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
cd $W
$E/python $W/scripts/fig2e_locus_check.py > $W/jobs/7020/fig2e.log 2>&1
$E/Rscript $W/scripts/fig2c_coding_owner.R > $W/jobs/7020/fig2c.log 2>&1
$E/Rscript $W/scripts/fig2b_same_row.R > $W/jobs/7020/fig2b.log 2>&1
$E/pip freeze > $W/pip_freeze.txt 2>/dev/null || true
