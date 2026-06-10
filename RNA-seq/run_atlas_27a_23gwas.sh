#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/27a_atlas_23gwas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/27a_atlas_23gwas_%j.out
# Phase 4a of the 23-GWAS COLOC refactor: rebuild the multi-evidence atlas
# so coloc_best_susie_pp4 / coloc_best_pp4 reflect the 23-study portfolio.
set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "host=$(hostname) start=$(date)"
Rscript 27a_assemble_evidence_atlas.R
echo "27a DONE rc=$? end=$(date)"
