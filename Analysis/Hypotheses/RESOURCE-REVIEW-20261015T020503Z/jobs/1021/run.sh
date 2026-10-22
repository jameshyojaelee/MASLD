#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1021
export TMPDIR=$J/tmp; mkdir -p $TMPDIR
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $J/fig2b_pip_variants.R
