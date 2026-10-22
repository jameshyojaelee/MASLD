#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4041
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp PYTHONDONTWRITEBYTECODE=1
mkdir -p $J/tmp $J/out
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$PY" -m pip freeze > $J/out/pip_freeze.txt 2>/dev/null || true
"$PY" $J/sim_pr.py $J/out 30605 5
echo DONE
