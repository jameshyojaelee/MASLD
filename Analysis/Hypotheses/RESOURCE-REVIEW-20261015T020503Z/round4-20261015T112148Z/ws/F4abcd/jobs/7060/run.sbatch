#!/bin/bash
#SBATCH --job-name=7060
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4abcd/jobs/7060/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4abcd/jobs/7060/%j.err
# R10-01: independent rerun of the review job 1101 Hotspot unweighted matched-random vocabulary null
# (arms A 469-record reproduction, B 364-unit, C 457-record Fig 4E roster). Script and lib copied unchanged.
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4abcd/jobs/7060
cd "$J"
export TMPDIR="$J/tmp" XDG_CACHE_HOME="$J/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript t3_hotspot_vocab_null_774.R
