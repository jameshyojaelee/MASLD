#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/2121
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp PYTHONNOUSERSITE=1
cd $J
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python refit.py 19999
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > pip_freeze.txt 2>/dev/null || true
