#!/bin/bash
#SBATCH --job-name=code-review --account=nslab --qos=nslab --partition=cpu,io
#SBATCH --cpus-per-task=1 --mem=4G --time=0:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/figures_candidate_20261102T160000Z/slurm_%j.out
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/figures_candidate_20261102T160000Z/make_panels.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/figures_candidate_20261102T160000Z
for f in /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/figures_candidate_20261102T160000Z/*.pdf; do pdftoppm -png -r 300 -singlefile $f /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/figures_candidate_20261102T160000Z/preview_$(basename $f .pdf); done
echo EXIT_OK
