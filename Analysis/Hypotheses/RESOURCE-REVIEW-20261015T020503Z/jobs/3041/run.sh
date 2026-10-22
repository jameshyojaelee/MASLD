#!/bin/bash
set -euo pipefail
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01
export TMPDIR=$S/tmp XDG_CACHE_HOME=$S/tmp MPLCONFIGDIR=$S/tmp PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3041/gwas_strength.py
