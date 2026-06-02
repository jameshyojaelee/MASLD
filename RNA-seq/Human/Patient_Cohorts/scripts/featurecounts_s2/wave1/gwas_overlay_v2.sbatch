#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=dream
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/gwas_overlay_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/gwas_overlay_v2_%j.err
# Wave1 [B] TERMINAL: 30_gwas_overlay_v2 (GWAS loci x DEG overlay). Reads canonical -s2 dream.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/gwas_spatial_convergence/logs
echo "=== gwas overlay v2 | $(date) ==="
Rscript 30_gwas_overlay_v2.R
if [ $? -ne 0 ]; then echo "FAILED: 30_gwas_overlay_v2"; exit 1; fi
echo "=== done $(date) ==="
