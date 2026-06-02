#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/fig1_compact_%j.out
#SBATCH --error=scripts/figures/logs/fig1_compact_%j.err
# Re-assemble fig1_compact.pdf only (after promoting the GWAS alluvial to panel d).
set +e
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE" || exit 1
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript scripts/figures/fig1_compact.R
echo "[exit fig1_compact=$?]"
ls -l --time-style=+%H:%M figures/main/fig1_atlas_overview/fig1_compact.pdf
