#!/bin/bash
#SBATCH --job-name=3921
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3921
mkdir -p $J/tmp
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
$PY $J/gate_c02_feasibility.py > $J/gate_c02_feasibility.log 2>&1
$PY -m pip freeze > $J/pip_freeze.txt 2>/dev/null || true
echo finished
