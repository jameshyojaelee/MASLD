#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4011
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $J/check.py
/gpfs/commons/home/jameslee/.local/bin/pdftoppm -r 300 -png -singlefile /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig7_model_framework/fig7E_verdict_vocabulary_and_refusal_thresholds.pdf $J/fig7E_render
ls -la $J
