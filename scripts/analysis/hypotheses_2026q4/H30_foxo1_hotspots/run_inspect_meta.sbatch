#!/bin/bash
#SBATCH --job-name=H30_meta
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/analysis/hypotheses_2026q4/H30_foxo1_hotspots/00_inspect_meta.R \
  RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds
