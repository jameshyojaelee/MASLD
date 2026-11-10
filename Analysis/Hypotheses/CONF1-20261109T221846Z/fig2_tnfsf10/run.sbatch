#!/bin/bash
#SBATCH --job-name=CONF1tnf
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/fig2_tnfsf10/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/fig2_tnfsf10/logs/%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1
STAGE=${1:-all}
if [ "$STAGE" = all ] || [ "$STAGE" = extract ]; then /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/fig2_tnfsf10/extract_tnfsf10.py; fi
if [ "$STAGE" = all ] || [ "$STAGE" = plot ]; then /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/fig2_tnfsf10/fig2_tnfsf10_locus.R; fi
