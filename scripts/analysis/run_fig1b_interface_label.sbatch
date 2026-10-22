#!/bin/bash
#SBATCH --job-name=fig1b_interface
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=logs/fig1b_interface_%j.out
#SBATCH --error=logs/fig1b_interface_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
FIG1B_INTERFACE_ROOT="figures/candidates/fig1b-native-interface_20261014_${SLURM_JOB_ID}"
export MASLD_PROJECT_ROOT="$PWD"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/figures/fig1b_current_interface.R "$FIG1B_INTERFACE_ROOT"
cp scripts/figures/fig1b_current_interface.R scripts/analysis/run_fig1b_interface_label.sbatch "$FIG1B_INTERFACE_ROOT/"
pdftotext "$FIG1B_INTERFACE_ROOT/fig1b_complementary_maps.pdf" "$FIG1B_INTERFACE_ROOT/printed_values.txt"
pdftoppm -r 150 -singlefile -png "$FIG1B_INTERFACE_ROOT/fig1b_complementary_maps.pdf" "$FIG1B_INTERFACE_ROOT/preview"
