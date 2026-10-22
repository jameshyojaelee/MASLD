#!/bin/bash
#SBATCH --job-name=2081
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/2081
mkdir -p $J/tmp
export TMPDIR=$J/tmp XDG_CACHE_HOME=$J/tmp MPLCONFIGDIR=$J/tmp OUT_DIR=$J
cd $J
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $J/rederive_fig2b.R > $J/rederive_fig2b.log 2>&1
