#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=fig1h
#SBATCH --output=logs/fig1h/fig1h_%j.log

set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript scripts/figures/fig1h_integration_discovery.R
