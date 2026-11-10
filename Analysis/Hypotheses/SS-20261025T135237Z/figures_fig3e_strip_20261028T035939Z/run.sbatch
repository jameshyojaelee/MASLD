#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/figures_fig3e_strip_20261028T035939Z/slurm_%j.out
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/figures_fig3e_strip_20261028T035939Z/make_strip.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/figures_fig3e_strip_20261028T035939Z
echo EXIT_OK
