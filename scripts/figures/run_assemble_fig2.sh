#!/bin/bash
#SBATCH --job-name=assemble_fig2
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/assemble_fig2_%j.out
#SBATCH --error=logs/assemble_fig2_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs figures/main/fig2_progression_sex

# Run the canonical May-17 compositor (do not modify the script itself).
Rscript scripts/figures/fig2_progression.R

# Canonical output is fig2_progression.pdf; mirror to fig2_composite.pdf
# (requested output name) so downstream consumers can find either.
cp -f figures/main/fig2_progression_sex/fig2_progression.pdf \
      figures/main/fig2_progression_sex/fig2_composite.pdf

ls -la figures/main/fig2_progression_sex/fig2_progression.pdf \
       figures/main/fig2_progression_sex/fig2_composite.pdf
echo "[done] fig2 composite assembled"
