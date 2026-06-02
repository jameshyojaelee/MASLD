#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=343i_fstage_ensemble
#SBATCH --output=Analysis/SingleCell/logs/343i_fstage_ensemble_%j.out
#SBATCH --error=Analysis/SingleCell/logs/343i_fstage_ensemble_%j.err
set -eo pipefail
mkdir -p Analysis/SingleCell/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rapids_singlecell python -u Analysis/SingleCell/scripts/343i_fstage_ensemble.py
# Plot via rnaseq R env
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/343i_fstage_ensemble_plot.R
