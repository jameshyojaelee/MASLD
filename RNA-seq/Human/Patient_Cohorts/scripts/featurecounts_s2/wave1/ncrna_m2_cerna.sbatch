#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/ncRNA_M2_cerna_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/ncRNA_M2_cerna_%j.err
# Wave1 [D] ncRNA M2: 54_cerna_network (after M1). Curated fallbacks + TargetScan cache.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "=== ncRNA M2 54_cerna_network | $(date) ==="
Rscript 54_cerna_network.R
if [ $? -ne 0 ]; then echo "FAILED: 54_cerna_network"; exit 1; fi
echo "=== done $(date) ==="
