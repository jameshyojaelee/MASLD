#!/bin/bash
#SBATCH --job-name=edgeR
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/edger_retest_%j.out
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/edger_retest_%j.err
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
# rnaseq env has edgeR + limma; ignore broken ~/.local numexpr-style shadowing.
export PYTHONNOUSERSITE=1
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
Rscript Analysis/ATAC/Human_Multiome/scripts/17_pseudobulk_de_retest.R
