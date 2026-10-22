#!/bin/bash
set -euo pipefail
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3131/tmp
mkdir -p "$TMPDIR"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 FLEXIBLAS_NUM_THREADS=4
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3131
GATE_B=${GATE_B:-400} /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript power.R
echo GATE_JOB_COMPLETE
