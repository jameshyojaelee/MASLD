#!/bin/bash
#SBATCH --job-name=figS05_atac_supp
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS05_atac_supp_%j.out
#SBATCH --error=logs/figS05_atac_supp_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/supplementary/figS05_epigenomic_spatial
Rscript scripts/figures/figS05_a_ld_aware_null_comparison.R
Rscript scripts/figures/figS05_b_motif_quality_tier.R
Rscript scripts/figures/figS05_c_gc_background_sensitivity.R
Rscript scripts/figures/figS05_d_stage_da_heatmap.R
