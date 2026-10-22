#!/bin/bash
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3821
export TMPDIR=$PWD/tmp; mkdir -p $TMPDIR
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/finemapping/bin/Rscript pilot_pnpla3.R
