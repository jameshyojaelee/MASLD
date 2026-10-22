#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/logs/%j.out
# AN3 task 2: Fig 3C four-cohort candidate (rerun after builder fix): source table, render, PNG proof.
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/code/00_build_fig3c_source_table.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/code/render_fig3c_four_cohort.R 2>&1 | tee /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/logs/render_fig3c.log
pdftoppm -r 300 -png -singlefile /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/panels/fig3c_four_cohort_overlap_bars.pdf /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/check/fig3c_four_cohort_overlap_bars
pdftotext -layout /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/panels/fig3c_four_cohort_overlap_bars.pdf /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/check/fig3c_four_cohort_overlap_bars.txt
env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3c-four-cohort-20261021T124839Z/logs/pip_freeze.txt
