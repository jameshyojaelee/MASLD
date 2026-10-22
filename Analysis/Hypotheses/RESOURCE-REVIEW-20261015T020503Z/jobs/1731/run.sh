#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1731
mkdir -p $J/tmp
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
$PY -m pip freeze > $J/out/pip_freeze.txt 2>/dev/null || true
$PY $J/rerun_spatial.py
