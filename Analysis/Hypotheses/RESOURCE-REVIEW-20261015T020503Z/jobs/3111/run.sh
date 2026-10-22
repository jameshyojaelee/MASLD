#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3111
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/he_molecular/bin/python
cd $J
$PY -m pip freeze > pip_freeze.txt 2>/dev/null || true
$PY morph.py 360
