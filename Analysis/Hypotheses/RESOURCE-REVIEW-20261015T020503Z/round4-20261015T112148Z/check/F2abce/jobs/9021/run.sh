#!/bin/bash
set -euo pipefail
CK=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2abce
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$CK/tmp
echo "start $(date -Is) on $(hostname)"
$ENV/python $CK/scripts/02_rederive.py $CK
$ENV/Rscript $CK/scripts/03_annotate.R $CK
$ENV/pip freeze > $CK/out/pip_freeze.txt 2>/dev/null || $ENV/python -m pip freeze > $CK/out/pip_freeze.txt
echo "done $(date -Is)"
