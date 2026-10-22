#!/bin/bash
# Gate L01-04 pilot: nucleus-level peak-gene association feasibility in GSE296875 hepatocytes.
set -euo pipefail
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-04
export TMPDIR=$S/tmp XDG_CACHE_HOME=$S/tmp MPLCONFIGDIR=$S/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
echo "host $(hostname) start $(date -Is)"
$PY -W ignore pilot_links.py
$PY -m pip freeze > pip_freeze.txt 2>/dev/null || true
echo "end $(date -Is)"
