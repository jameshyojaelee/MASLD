#!/bin/bash
set -euo pipefail
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1111/tmp
mkdir -p $TMPDIR
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript perm_limma.R
