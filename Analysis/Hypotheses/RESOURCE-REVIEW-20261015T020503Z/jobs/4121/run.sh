#!/usr/bin/env bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4121
T=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify2-R14-02/tmp
mkdir -p "$T"
export TMPDIR="$T" XDG_CACHE_HOME="$T" MPLCONFIGDIR="$T"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$J"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -B trait_split.py > trait_split.log 2>&1
echo "exit $?"
