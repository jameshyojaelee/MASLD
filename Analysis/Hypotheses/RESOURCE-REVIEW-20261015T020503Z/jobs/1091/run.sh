#!/bin/bash
set -euo pipefail
export PYTHONNOUSERSITE=1
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1091/tmp
mkdir -p "$TMPDIR"
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1091
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python fig7b_mpra_multiplicity.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > pip_freeze.txt 2>/dev/null || true
