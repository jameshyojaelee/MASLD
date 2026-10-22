#!/bin/bash
#SBATCH --job-name=7520
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/jobs/7520/%j.out
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $D/scripts/fig2b_multirow_check.R
