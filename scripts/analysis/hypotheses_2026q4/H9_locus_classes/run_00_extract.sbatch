#!/bin/bash
#SBATCH --job-name=H9_extract
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H9-20261009T140001Z/logs/extract-%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/analysis/hypotheses_2026q4/H9_locus_classes/00_extract_supplement.R Analysis/Hypotheses/H9-20261009T140001Z/downloads_untrusted/jamialahmadi2024_nat_med/jamialahmadi2024_MOESM3.xlsx Analysis/Hypotheses/H9-20261009T140001Z/inputs/jamialahmadi2024_tsv
