#!/bin/bash
#SBATCH --job-name=7560
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F4abcd/jobs/7560/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F4abcd/jobs/7560/%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/producers
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript 02r_render_fig4c_h59.R > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/producers/02r_render.log 2>&1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair
pdftotext -layout fig4c_hotspot_stage_heatmap_H59_candidate_repair.pdf png_check/fig4c_repair.txt
pdftotext -layout figs4c_hotspot_stage_heatmap_full_H59_candidate_repair.pdf png_check/figs4c_repair.txt
pdftoppm -png -r 200 fig4c_hotspot_stage_heatmap_H59_candidate_repair.pdf png_check/fig4c_repair
pdftoppm -png -r 200 figs4c_hotspot_stage_heatmap_full_H59_candidate_repair.pdf png_check/figs4c_repair
