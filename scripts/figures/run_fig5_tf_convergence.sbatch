#!/bin/bash
#SBATCH --job-name=fig5_tf_conv
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig5_tf_convergence_%j.out
#SBATCH --error=logs/fig5_tf_convergence_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig5_convergence/panels

Rscript scripts/figures/fig5_panel_tf_convergence_atac_rna_coloc.R
