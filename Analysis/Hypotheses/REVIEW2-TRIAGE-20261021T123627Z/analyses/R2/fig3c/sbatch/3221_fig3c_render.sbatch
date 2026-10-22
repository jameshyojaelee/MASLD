#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/R2/fig3c/logs/%j.out
# R2 fig3c: Fig 3C five-cohort (GSE126848 filtered refit) candidate: source table, render, PNG proof, pdftotext.
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/code/00_build_fig3c_source_table.py | tee /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/logs/build_fig3c.log
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/code/render_fig3c_five_filtered.R 2>&1 | tee /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/logs/render_fig3c.log
pdftoppm -r 300 -png -singlefile /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/panels/fig3c_cohort_overlap_bars.pdf /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/check/fig3c_cohort_overlap_bars
pdftotext -layout /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/panels/fig3c_cohort_overlap_bars.pdf /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/check/fig3c_cohort_overlap_bars.txt
env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-five-filtered-20261021T172748Z/logs/pip_freeze.txt
