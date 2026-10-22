#!/bin/bash
#SBATCH --job-name=7095
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/logs/7095_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/logs/7095_%j.err
set -euo pipefail
WS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef
CAND=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F5ef
export TMPDIR=$WS/tmp XDG_CACHE_HOME=$WS/tmp R_LIBS_USER=
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$ENV/Rscript $WS/scripts/04_render_fig5e_candidate.R
pdftotext -layout $CAND/fig5e_program_context_117_candidate.pdf $WS/out/fig5e_candidate_pdftotext.txt
pdftoppm -png -r 300 -singlefile $CAND/fig5e_program_context_117_candidate.pdf $CAND/fig5e_program_context_117_candidate_check
