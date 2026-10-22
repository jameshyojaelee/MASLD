#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1041
export PYTHONNOUSERSITE=1 TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY $J/perm_and_content.py $J/out 9999 200
$PY -m pip freeze > $J/out/pip_freeze.txt 2>/dev/null || true
