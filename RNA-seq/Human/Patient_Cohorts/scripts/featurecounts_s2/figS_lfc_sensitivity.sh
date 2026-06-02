#!/bin/bash
#SBATCH --job-name=fgsea
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_figS_lfc_sensitivity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_figS_lfc_sensitivity_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== figS_lfc_sensitivity (atlas Hallmark ORA, -s 2) ==="
echo "Started: $(date)"

Rscript scripts/figures/figS_lfc_sensitivity.R
if [ $? -ne 0 ]; then echo "FAILED"; exit 1; fi

echo "Finished: $(date)"
