#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/zonation_43_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/zonation_43_%j.err
# Wave1 [G] ATLAS FEEDER: 43_zonation_classification -> deg_zonation_classification.csv (27a:1907).
# Reads current -s2 multi_evidence_atlas.csv for dream_logFC/padj columns (atlas already -s2, May 28).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 43 zonation classification | $(date) ==="
Rscript RNA-seq/43_zonation_classification.R
if [ $? -ne 0 ]; then echo "FAILED: 43_zonation_classification"; exit 1; fi
echo "=== done $(date) ==="
