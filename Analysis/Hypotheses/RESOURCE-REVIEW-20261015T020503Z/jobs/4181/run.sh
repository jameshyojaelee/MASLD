#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4181
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -m pip freeze > $J/pip_freeze.txt 2>/dev/null || true
/usr/bin/time -v $PY $J/verify_g1_01.py $J > $J/verify_g1_01.log 2> $J/verify_g1_01.time
